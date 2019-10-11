/* Copyright (C) 2018-2019 Bho Matthiesen, Christoph Hellings
 * 
 * This program is used in the article:
 *
 * Bho Matthiesen, Christoph Hellings, Eduard A. Jorswieck, and Wolfgang
 * Utschick, "Mixed Monotonic Programming for Fast Global Optimization,"
 * submitted to IEEE  Transactions on Signal Processing.
 * 
 * 
 * License:
 * This program is licensed under the GPLv2 license. If you in any way use this
 * code for research that results in publications, please cite our original
 * article listed above.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. */


#include <iostream>
#include <string>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <sstream>

#include "PA.h"
#include "H5Cpp.h"

template <size_t Dim, size_t outDim> class result;

template <size_t Dim>
struct input
{
	size_t wpidx;
	size_t cidx;
	size_t pidx;

	double ck[Dim];
	bool beta[Dim][Dim];
	double Rmin[Dim];

	input(const char *name, const size_t ppidx);

private:
	template <size_t, size_t> friend class result;
	input() { }; // creates uninitalized
};

template <size_t Dim>
input<Dim>::input(const char *name, const size_t wpidx)
	: wpidx(wpidx)
{
	H5::H5File h5(name, H5F_ACC_RDONLY);

	// alpha
	{
		auto ds = h5.openDataSet("ck");
		auto sp = ds.getSpace();

		const hsize_t count[] = {1, Dim};
		const hsize_t start[] = {wpidx, 0, 0};
		sp.selectHyperslab(H5S_SELECT_SET, count, start);

		H5::DataSpace memspace(1, count+1);
		memspace.selectAll();

		ds.read(&ck, H5::PredType::NATIVE_DOUBLE, memspace, sp);
	}

	// beta
	{
		auto ds = h5.openDataSet("beta");
		auto sp = ds.getSpace();

		const hsize_t count[] = {1, Dim, Dim};
		const hsize_t start[] = {wpidx, 0, 0};
		sp.selectHyperslab(H5S_SELECT_SET, count, start);

		H5::DataSpace memspace(2, count+1);
		memspace.selectAll();

		ds.read(&beta, H5::PredType::NATIVE_HBOOL, memspace, sp);
	}

	// Rmin
	{
		auto ds = h5.openDataSet("Rmin");
		auto sp = ds.getSpace();

		const hsize_t count[] = {1, Dim};
		const hsize_t start[] = {wpidx, 0};
		sp.selectHyperslab(H5S_SELECT_SET, count, start);

		H5::DataSpace memspace(1, count+1);
		memspace.selectAll();

		ds.read(&Rmin, H5::PredType::NATIVE_DOUBLE, memspace, sp);
	}
}

template <size_t Dim>
std::ostream& operator<< (std::ostream &out, const input<Dim>& data)
{
	out << "ck: [ ";
	for (size_t i = 0; i < Dim; ++i)
		out << data.ck[i] << " ";
	out << "]" << std::endl;

	out << "Rmin: [ ";
	for (size_t i = 0; i < Dim; ++i)
		out << data.Rmin[i] << " ";
	out << "]" << std::endl;

	out << "beta:  [[ ";
	for (size_t i = 0; i < Dim; ++i)
	{
		if (i > 0)
			out << std::endl << "        [ ";

		for (size_t j = 0; j < Dim; ++j)
			out << data.beta[i][j] << " ";
		out << "]";
	}
	out << "]" << std::endl;

	return out;
}


template <size_t Dim, size_t outDim = Dim>
class result
{
public:
	struct RT
	{
		double objval;
		double xopt[outDim];
		std::string status;
		unsigned long long iterations;
		unsigned long long lastUpdate;
		double epsilon;
		bool useRelTol;
		double runtime;
		unsigned long long maxQueueSize;
		unsigned dataSize;

		template <class T> RT(const T& o);

	private:
		friend class result<Dim,outDim>;
		RT() {};
	};

	// ctor
	result(const std::string& name, const input<Dim>& data);
	result(const result& ) = delete;
	result& operator=(const result& ) = delete;

	void writeResult(const std::string& name, const RT& res);
	std::unique_ptr<RT> getResult(const std::string& name) const;

private:
	std::unique_ptr<H5::H5File> h5;
	H5::CompType mtRT;
};


template <size_t Dim, size_t outDim>
template <class T>
result<Dim, outDim>::RT::RT(const T& o)
{
	static_assert(o.xopt.size() == outDim);

	objval = o.optval;
	std::copy(o.xopt.begin(), o.xopt.end(), xopt);
	status = o.statusStr;
	iterations = o.iter;
	lastUpdate = o.lastUpdate;
	epsilon = o.getEpsilon();

	if constexpr(std::is_base_of<pa::PA<outDim>,T>::value)
		useRelTol = false;
	else
		useRelTol = o.useRelTol;

	runtime = o.runtime;
	maxQueueSize = o.max_queue_size;
	dataSize = o.data_size;
}


template <size_t Dim, size_t outDim>
result<Dim,outDim>::result(const std::string& name, const input<Dim>& data)
	: mtRT(sizeof(RT))
{
	// open file
	bool init = false;

	try
	{
		h5 = std::make_unique<H5::H5File>(name, H5F_ACC_EXCL);
		init = true;
	}
	catch (H5::FileIException& e)
	{
		try
		{
			h5 = std::make_unique<H5::H5File>(name, H5F_ACC_RDWR);
			init = false;
		}
		catch (H5::FileIException& e)
		{
			h5 = std::make_unique<H5::H5File>(name, H5F_ACC_TRUNC);
			init = true;
		}
	}

	// create memory datatype for input<>
	H5::CompType mtInput(sizeof(input<Dim>));
	mtInput.insertMember("wpidx", HOFFSET(input<Dim>, wpidx), H5::PredType::NATIVE_HSIZE);
	{
		const hsize_t dims[] = {Dim};
		H5::ArrayType dta(H5::PredType::NATIVE_DOUBLE, 1, dims);
		mtInput.insertMember("ck", HOFFSET(input<Dim>, ck), dta);
		mtInput.insertMember("Rmin", HOFFSET(input<Dim>, Rmin), dta);
	}
	{
		const hsize_t dims[] = {Dim, Dim};
		H5::ArrayType dta(H5::PredType::NATIVE_HBOOL, 2, dims);
		mtInput.insertMember("beta", HOFFSET(input<Dim>, beta), dta);
	}

	// initalize memory datatype for RT
	mtRT.insertMember("objval", HOFFSET(RT, objval), H5::PredType::NATIVE_DOUBLE);
	{
		const hsize_t dims[] = {outDim};
		H5::ArrayType dta(H5::PredType::NATIVE_DOUBLE, 1, dims);
		mtRT.insertMember("xopt", HOFFSET(RT, xopt), dta);
	}
	mtRT.insertMember("status", HOFFSET(RT, status), H5::StrType(0, H5T_VARIABLE));
	mtRT.insertMember("iterations", HOFFSET(RT, iterations), H5::PredType::NATIVE_ULLONG);
	mtRT.insertMember("lastUpdate", HOFFSET(RT, lastUpdate), H5::PredType::NATIVE_ULLONG);
	mtRT.insertMember("epsilon", HOFFSET(RT, epsilon), H5::PredType::NATIVE_DOUBLE);
	mtRT.insertMember("useRelTol", HOFFSET(RT, useRelTol), H5::PredType::NATIVE_HBOOL);
	mtRT.insertMember("runtime", HOFFSET(RT, runtime), H5::PredType::NATIVE_DOUBLE);
	mtRT.insertMember("maxQueueSize", HOFFSET(RT, maxQueueSize), H5::PredType::NATIVE_ULLONG);
	mtRT.insertMember("dataSize", HOFFSET(RT, dataSize), H5::PredType::NATIVE_UINT);

	if (init) // initalize
	{
		auto ds = h5->createDataSet("input", mtInput, {});
		ds.write(&data, mtInput);
	}
	else // validate file
	{
		input<Dim> tmp;
		h5->openDataSet("input").read(&tmp, mtInput);

		if (tmp.wpidx != data.wpidx)
			throw std::logic_error("wpidx mismatch in existing result file");
	}
}


template <size_t Dim, size_t outDim>
void
result<Dim, outDim>::writeResult(const std::string& name, const RT& res)
{
	auto ds = h5->createDataSet(name, mtRT, {});
	ds.write(&res, mtRT);
}


template <size_t Dim, size_t outDim>
std::unique_ptr<typename result<Dim,outDim>::RT>
result<Dim,outDim>::getResult(const std::string& name) const
{
	std::unique_ptr<RT> res(new RT()); // make_unique not applicable due to private ctor

	auto ds = h5->openDataSet(name);
	ds.read(res.get(), mtRT);

	return res;
}


template <template <size_t> class T, size_t Dim, size_t outDim = Dim>
void
aloha_helper(const char* infn, const unsigned wpidx, const char* name, bool reduction = false)
{
	// get input data
	//H5::Exception::dontPrint();
	const input<Dim> data(infn, wpidx);

	// create output file
	std::stringstream ss;

	const char *sd = std::getenv("JOB_HPC_SAVEDIR");
	if (sd != nullptr)
	{
		ss << sd;

		if (sd[std::strlen(sd)-1] != '/')
			ss << "/";
	}

	ss << "result_" << name << "_" << wpidx << ".h5";
	std::cout << ss.str() << std::endl;
	result<Dim,outDim> resF(ss.str(), data);

	// set up problem
	T<Dim> aloha;

	aloha.setPrecision(1e-3);
	aloha.useRelTol = false;
	aloha.disableReduction = !reduction;
	aloha.outputEvery = 1e6;
	//aloha.enablePruning = false;

	for (size_t i = 0; i < Dim; ++i)
	{
		aloha.ck[i] = data.ck[i];
		aloha.Rmin[i] = data.Rmin[i];

		for (size_t j = 0; j < Dim; ++j)
			aloha.beta[i][j] = data.beta[i][j];
	}

	std::cout << data << std::endl;

	// optimize
	aloha.optimize();
	aloha.optval = std::exp(aloha.optval);

	// save
	typename decltype(resF)::RT res(aloha);
	resF.writeResult(name, res);

	// debug output
	for (size_t i = 0; i < Dim; ++i) {
		double tmp = 1.0;
		for (size_t j = 0; j < Dim; ++j)
			if (aloha.beta[i][j] > 0)
				tmp *= (1.0 - aloha.xopt[j]);

		std::cout << "R_" << i << " = " << (aloha.ck[i] * aloha.xopt[i] * tmp) << ", ";
	}
	std::cout << std::endl;
}


bool
stob(const std::string& s)
{
	int tmp = std::stoi(s);
	switch (tmp)
	{
		case 0:
			return false;
		case 1:
			return true;
		default:
			throw std::invalid_argument("stob");
	}
}
