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


#ifndef _ALOHA_GP_H_
#define _ALOHA_GP_H_

#include <array>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <sstream>

extern "C" {
	#include "mosek.h"
	#include "dgopt.h"
}

#include "Mosek.h"
#include "PA.h"
#include "aloha_helper.h"


using namespace pa;

template <size_t Dim>
class ALOHA : public PA<Dim>
{
	using typename PA<Dim>::vtype;
	
	public:
		ALOHA() : PA<Dim>() {}

	private:
		vtype projection(const vtype& x) const override;
};


template <size_t Dim>
struct userdata
{
	std::array<double, Dim> Rmin;
	std::array<double, Dim> ck;

	// Mosek
	MSKtask_t task = nullptr;
	double v[3*Dim]; // numVars
	dgohand_t nlh = nullptr;

	constexpr size_t dim() { return Dim; }

	userdata(const bool b[][Dim]);

	~userdata() { destroy(); };

private:
	void destroy();
};

template <size_t Dim>
void
userdata<Dim>::destroy()
{
	if (task)
	{
		MSK_freedgo(task, &nlh);
		MSK_deletetask(&task);
	}
}

template <size_t Dim>
userdata<Dim>::userdata(const bool b[][Dim]) try
{
	const size_t numVars = 3*Dim;
	const size_t numConstr = 2*Dim+1; // A*delta=0
	// init Mosek
	mskcall(MSK_maketask(MosekEnv::getEnv(), numConstr, numVars, &task));

	mskcall(MSK_appendvars(task, numVars));
	mskcall(MSK_appendcons(task, numConstr));
	mskcall(MSK_putintparam(task, MSK_IPAR_INTPNT_MULTI_THREAD, MSK_OFF));


	/* set all variables positive */
	for (size_t i = 0; i < numVars; ++i)
		mskcall(MSK_putvarbound(task, i, MSK_BK_LO, 0.0, +MSK_INFINITY));

	/* A "entries" for obj == bound */
	mskcall(MSK_putconbound(task, 0, MSK_BK_FX, 1.0, 1.0));
	for (size_t i = 1; i < numConstr; ++i)
		mskcall(MSK_putconbound(task, i, MSK_BK_FX, 0.0, 0.0));

	/* A entries for (7.8b) */
	for (size_t i = 0; i < Dim; ++i)
	{
		mskcall(MSK_putaij(task, 0, i, 1.0)); // alpha
		mskcall(MSK_putaij(task, 1+i, i, -1.0)); // theta_i

		for (size_t j = 0; j < Dim; ++j)
		{
			if (b[i][j])
				mskcall(MSK_putaij(task, 1+Dim+j, i, -1.0)); // \hat theta_i
		}
	}

	/* A entries for (7.8c) */
	for (size_t i = 0; i < Dim; ++i)
	{
		v[Dim+2*i] = 1.0;
		v[Dim+2*i+1] = 1.0;

		mskcall(MSK_putaij(task, 1+i, Dim + 2*i, 1.0)); // theta_i
		mskcall(MSK_putaij(task, 1+Dim+i, Dim + 2*i+1, 1.0)); // \hat theta_i
	}
	
	/* init p */
	MSKint32t p[2*Dim+1];
	p[0] = 0;
	for (size_t i = 0; i < Dim; ++i)
	{
		p[i+1] = 1;
		p[i+Dim+1] = 2;
	}

	mskcall(MSK_dgosetup(task, numVars, numConstr , 2*Dim+1, v, p, &nlh));
}
catch (...)
{
	destroy();
	throw;
}

template <size_t Dim>
typename ALOHA<Dim>::vtype
ALOHA<Dim>::projection(const typename ALOHA<Dim>::vtype& x) const
{
	userdata<Dim>* ud = static_cast<userdata<Dim>*>(this->user_data);

	/* update v entries for (7.8b) */
	for (size_t i = 0; i < Dim; ++i)
		ud->v[i] = x[i]/ud->ck[i];

	mskcall(MSK_dgoupdatev(&ud->nlh, ud->v));
	mskcall(MSK_optimize(ud->task));

	/* obtain solution */
	double alpha;
	MSK_getsucslice(ud->task, MSK_SOL_ITR, 0, 1, &alpha);
	alpha = std::exp(alpha);

	vtype ret;
	for (size_t i = 0; i < Dim; ++i)
		ret[i] = x[i] * alpha;

	return ret;
}


template <size_t Dim>
bool
inH(const typename pa::PA<Dim>::vtype& x, void* user_data)
{
	userdata<Dim>* ud = (userdata<Dim> *) user_data;
	
	bool ret = true;

	for (size_t i = 0; i < Dim; ++i)
		ret = ret && (x[i] >= ud->Rmin[i]);

	return ret;
}


template <size_t Dim>
double
obj(const typename pa::PA<Dim>::vtype& x, void* user_data)
{
	double ret = 1; // prop fair utility

	for (size_t i = 0; i < Dim; ++i)
		ret *= (x[i]);  // prop fair utility

	return std::log(ret);
}


template <size_t Dim>
bool
inG(const typename pa::PA<Dim>::vtype& x, void* user_data)
{
	userdata<Dim>* ud = (userdata<Dim> *) user_data;
	
	/* A entries for (7.8b) */
	for (size_t i = 0; i < Dim; ++i)
		ud->v[i] = x[i]/ud->ck[i];

	mskcall(MSK_dgoupdatev(&ud->nlh, ud->v));
	mskcall(MSK_optimize(ud->task));

	/* obtain solution */
	double alpha;
	MSK_getsucslice(ud->task, MSK_SOL_ITR, 0, 1, &alpha);

	return alpha > 0.0;
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
	result<Dim,Dim> resF(ss.str(), data);

	// set up problem
	T<Dim> aloha;
	userdata<Dim> ud(data.beta);

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
		aloha.setUB(i, data.ck[i]);
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

		std::cout << "R_" << i << " = " << aloha.xopt[i] << ", ";
	}
	std::cout << std::endl;
}
}

#endif
