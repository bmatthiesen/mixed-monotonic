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


#include <numeric>
#include <cstring>

#include <iostream>
#include <sstream>

extern "C" {
	#include "mkl.h"
}

#include "PA.h"
#include "aloha_helper.h"

using namespace pa;

constexpr size_t DIM = 6;


template <size_t Dim>
class ALOHA : public PA<2*Dim>
{
	public:
		using typename PA<2*Dim>::vtype;
		ALOHA() : PA<2*Dim>() {}

	private:
		vtype projection(const vtype& x) const override;
};

template <size_t Dim>
typename ALOHA<Dim>::vtype
ALOHA<Dim>::projection(const typename ALOHA<Dim>::vtype& x) const
{
	double beta = 1.0;
	vtype ret;

	for (size_t i = 0; i < Dim; ++i)
	{
		beta = std::min(beta, 1/(x[i] + x[Dim+i]));
	}
	for (size_t i = 0; i < 2*Dim; ++i)
	{
		ret[i] = beta*x[i];
	}

	return ret;
}





template <size_t Dim>
struct userdata
{
	double Rmin[Dim];
	double ck[Dim];
	bool beta[Dim][Dim];

	constexpr size_t dim() { return Dim; }
};


template <size_t Dim>
bool
inH(const typename ALOHA<Dim>::vtype& x, void* user_data)
{
	userdata<Dim>* ud = (userdata<Dim> *) user_data;
	
	/* check that all are >= 0.0 */
	bool ret = true;
	for (auto it = x.begin(); ret && it != x.end(); ++it)
		ret = ret && *it >= 0.0;

	for (size_t i = 0; i < Dim; ++i)
	{
		double tmp = 1.0;
		for (size_t j = 0; j < Dim; ++j)
			if (ud->beta[i][j] > 0)
				tmp *= x[Dim+j];
				
		tmp *= ud->ck[i] * x[i];
		if (tmp < ud->Rmin[i])
		{
			ret = false;
			break;
		}
	}
	
	return ret;
}


template <size_t Dim>
double
obj(const typename ALOHA<Dim>::vtype& x, void* user_data)
{
	userdata<Dim>* ud = (userdata<Dim> *) user_data;

	double ret = 1; // prop fair utility

	for (size_t i = 0; i < Dim; ++i) {
		double tmp = 1.0;
		for (size_t j = 0; j < Dim; ++j)
			if (ud->beta[i][j])
				tmp *= x[Dim+j];
		
		ret *= (ud->ck[i] * x[i] * tmp);  // prop fair utility
	}

	return std::log(ret);
}


template <size_t Dim>
bool
inG(const typename ALOHA<Dim>::vtype& x, void*)
{
	bool ret = true;
	for (size_t i = 0; i < Dim; ++i)
	{
		bool tmp = x[i] + x[Dim+i] <= 1;

		if (!tmp)
		{
			ret = false;
			break;
		}
	}

	return ret;
}

namespace pa {
	template <template <size_t> class T, size_t Dim>
	void
	aloha_helper(const char* infn, const unsigned wpidx, const char* name)
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
		result<Dim,2*Dim> resF(ss.str(), data);

		// set up problem
		T<Dim> aloha;
		userdata<Dim> ud;

		aloha.setPrecision(1e-3);
		aloha.setObjective(obj<Dim>);
		aloha.setInH(inH<Dim>);
		aloha.setInG(inG<Dim>);
		aloha.setUserData(&ud);
		//aloha.setShift(); // shift not necessary for this problem
		aloha.output_every = 1e2;

		for (size_t i = 0; i < Dim; ++i)
		{
			ud.ck[i] = data.ck[i];
			ud.Rmin[i] = data.Rmin[i];
			aloha.setUB(i, 1);
			aloha.setUB(i+Dim, 1);

			for (size_t j = 0; j < Dim; ++j)
				ud.beta[i][j] = data.beta[i][j];
		}

		std::cout << data << std::endl;

		// optimize
		aloha.optimize();
		aloha.optval = std::exp(aloha.optval);

		// save
		typename decltype(resF)::RT res(aloha);
		resF.writeResult(name, res);

		// debug output
		for (size_t i = 0; i < Dim; ++i)
		{
			double tmp = 1.0;

			for (size_t j = 0; j < Dim; ++j)
				if (ud.beta[i][j])
					tmp *= (1.0 - aloha.xopt[j]);

			std::cout << "R_" << i << " = " << (ud.ck[i] * aloha.xopt[i] * tmp) << ", ";
		}
		std::cout << std::endl;
	}
}


int
main(int argc, char *argv[]) try
{
	const int Dim = 3;

	// parse input arguments
	if (argc != 3)
		throw std::invalid_argument("");

	unsigned wpidx = std::stoi(argv[2]);

	pa::aloha_helper<ALOHA,Dim>(argv[1], wpidx, "PA");

	return 0;
}
catch (H5::Exception& e)
{
	std::cerr << "HDF Exception" << std::endl;
	e.printErrorStack();
	return -1;
}
catch(std::invalid_argument& e)
{
	std::cerr << "Usage: " << argv[0] << " wpfile wpidx" << std::endl;

	return -2;
}
