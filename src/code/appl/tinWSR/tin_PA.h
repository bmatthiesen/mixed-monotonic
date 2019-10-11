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

#ifndef _TIN_PA_H_
#define _TIN_PA_H_

#include <limits>
#include <array>
#include <numeric>
#include <cmath>

#include "PA.h"

template <size_t Dim>
struct userdata
{
	double N;
	std::array<double, Dim> Pmax;
	std::array<double, Dim> sigma;
	std::array<double, Dim> alpha;
	std::array<double, Dim> beta[Dim];

	double tmax;
	double tmin;

	constexpr size_t dim() { return Dim; }

	void preprocess();
};

template <size_t Dim>
bool
inH(const std::array<double, Dim>& x, void *user_data)
{
	userdata<Dim-1>* ud = (userdata<Dim-1> *) user_data;

	/* check that all are >= 0.0 */
	for (size_t i = 0; i < Dim-1; ++i)
		if (x[i] < 0.0)
			return false;

	return x[Dim-1] >= ud->tmin;
}

template <size_t Dim>
double
obj_f(const std::array<double, Dim>& x, void* user_data)
{
	userdata<Dim-1>* ud = (userdata<Dim-1> *) user_data;

	double ret = 1;
	
	for (size_t i = 0; i < ud->dim(); ++i)
	{
		double tmp = ud->alpha[i] * x[i] + std::inner_product(ud->beta[i].begin(), ud->beta[i].end(), x.begin(), ud->sigma[i]);

		if (tmp < 0)
			return -std::numeric_limits<double>::infinity();

		ret *= tmp;
	}

	return std::log2(ret);
}

template <size_t Dim>
double
obj_g(const std::array<double, Dim>& x, void* user_data)
{
	userdata<Dim-1>* ud = (userdata<Dim-1> *) user_data;

	double ret = 1;
	
	for (size_t i = 0; i < ud->dim(); ++i)
	{
		double tmp = std::inner_product(ud->beta[i].begin(), ud->beta[i].end(), x.begin(), ud->sigma[i]);

		if (tmp < 0)
			return -std::numeric_limits<double>::infinity();

		ret *= tmp;
	}

	return std::log2(ret);
}

template <size_t Dim>
bool
inG(const std::array<double, Dim>& x, void* user_data)
{
	userdata<Dim-1>* ud = (userdata<Dim-1> *) user_data;

	for (size_t i = 0; i < ud->dim(); ++i)
		if (x[i] > ud->Pmax[i])
			return false;

	return x[Dim-1] <= ud->tmax && x[Dim-1] + obj_g<Dim>(x, user_data) <= 0;
}

template <size_t Dim>
void
userdata<Dim>::preprocess()
{
	std::array<double, Dim+1> zero;
	zero.fill(0);

	std::array<double, Dim+1> max;
	for (size_t i = 0; i < Dim; ++i)
		max[i] = Pmax[i];
	max[Dim] = 0;

	tmin = -obj_g<Dim+1>(max, this);
	tmax = -obj_g<Dim+1>(zero, this);
}

template <size_t Dim>
double
obj(const std::array<double, Dim>& x, void* user_data)
{
	return x[Dim-1] + obj_f<Dim>(x, user_data);
}

template <size_t Dim>
struct TIN : public pa::PA<Dim+1>
{
	TIN() : pa::PA<Dim+1>()
	{
		this->setPrecision(1e-2);
		this->setObjective(obj<Dim+1>);
		this->setInH(inH<Dim+1>);
		this->setInG(inG<Dim+1>);
	}
};

#endif
