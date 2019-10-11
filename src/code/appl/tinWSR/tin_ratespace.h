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


#ifndef _TIN_RATESPACE_H
#define _TIN_RATESPACE_H

#include <array>
#include <numeric>

extern "C" {
	#include "mkl.h"
}

#include "PA.h"

template <size_t Dim>
struct userdata
{
	double N;
	std::array<double, Dim> Pmax;
	std::array<double, Dim> sigma;
	std::array<double, Dim> alpha;
	std::array<std::array<double, Dim>, Dim> beta;

	// computed by preprocess()
	std::array<std::array<double, Dim>, Dim> G;
	std::array<double, Dim> eta;

	constexpr size_t dim() { return Dim; }
	void preprocess(void); // call after setting alpha, beta and sigma
};

template <size_t Dim>
void
userdata<Dim>::preprocess(void)
{
	for (size_t i = 0; i < Dim; ++i)
	{
		for (size_t j = 0; j < Dim; ++j)
		{
			if (i == j)
			{
				G[i][j] = 0.0;

				if (beta[i][i] != 0.0)
					throw std::domain_error("beta_ii must be zero");
			}
			else
				G[i][j] = beta[i][j] / alpha[i];
		}

		eta[i] = sigma[i] / alpha[i];
	}
}

template <size_t Dim>
bool
inH(const std::array<double, Dim>& x, void*)
{
	/* check that all are >= 0.0 */
	bool ret = true;
	for (const auto e : x)
		ret = ret && e >= 0.0;
	
	return ret;
}

template <size_t Dim, bool propFair = false>
double
obj(const std::array<double, Dim>& x, void*)
{
	if constexpr(propFair)
	{
		double ret = 1;
		for (size_t i = 0; i < Dim; ++i)
		{
			ret *= x[i];
		}
		return ret;
	}
	else
		return std::accumulate(x.begin(), x.end(), 0.0);
}

template <size_t Dim>
bool
inG(const std::array<double, Dim>& x, void* user_data)
{
	userdata<Dim>* ud = (userdata<Dim> *) user_data;

	double G[Dim][Dim] __attribute__((aligned(64)));
	double eta[Dim] __attribute__((aligned(64)));

	for (size_t i = 0; i < Dim; ++i)
	{
		double gamma = std::pow(2.0, std::max(x[i],0.0)) - 1.0;

		eta[i] = gamma * ud->eta[i];

		for (size_t j = 0; j < Dim; ++j)
		{
			if (i == j)
				G[i][i] = 1.0;
			else
				G[i][j] = -gamma * ud->G[i][j];
		}
	}

	{
		lapack_int ipiv[Dim] __attribute__((aligned(64)));

		// WARNING: eta is overwritten, G is overwritten if info != 0
		auto info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, Dim, 1, &G[0][0], Dim, ipiv, eta, 1);

		if (info != 0)
			throw std::runtime_error("linear system solver failed with status " + std::to_string(info));
	}

	for (auto&& e : eta)
	{
		e = (abs(e)<1e-14) ? 0 : e;
	}

	bool ret = true;
	for (size_t i = 0; i < Dim; ++i)
	{
		bool tmp = eta[i] >= 0 && eta[i] <= ud->Pmax[i];

		if (!tmp)
		{
			ret = false;
			break;
		}
	}

	return ret;

	/*
	gamma = 2**(x) - 1;
	G_{kj} = gamma_k * beta_{kj} / alpha{k};
	G_{jj} = 0;

	eta_k = gamma_k * sigma[k] / alpha[k];

	Solve (I-G) *p = eta for p
	
	return all(0 <= p <= pmax)
	*/
}

template <size_t Dim, bool propFair = false>
struct Tin_Ratespace : public pa::PA<Dim>
{
	Tin_Ratespace(void *ud);
};

template <size_t Dim, bool propFair>
Tin_Ratespace<Dim,propFair>::Tin_Ratespace(void *ud)
 : pa::PA<Dim>()
{
	this->setPrecision(1e-2);
	this->setObjective(obj<Dim,propFair>);
	this->setInH(inH<Dim>);
	this->setInG(inG<Dim>);
	this->setShift();
	this->setUserData(ud);
}

#endif
