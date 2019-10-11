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


#ifndef _ALOHA_H_
#define _ALOHA_H_

#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>

#include "MMP.h"

using std::array;

template <size_t Dim>
class ALOHA : public MMPconstraints<Dim>
{
	using typename MMPconstraints<Dim>::vtypeS;
	using typename MMPconstraints<Dim>::RBox;

	public:
		double ck[Dim];
		bool beta[Dim][Dim];
		double Rmin[Dim];
		
		ALOHA();

		bool isOptimal() { return this->status == BRB<Dim>::Status::Optimal; }; // helper for python

	private:
		double MMPobj(const vtypeS& x, const vtypeS& y) const override;
		bool constraints(const vtypeS& x, const vtypeS& y) const override;
		vtypeS feasiblePoint(const RBox& r) const override
			{ return r.lb(); }
};


template <size_t D>
ALOHA<D>::ALOHA()
	: MMPconstraints<D>()
{
	this->setUB(1);
	this->setLB(0);
}


template <size_t D>
double
ALOHA<D>::MMPobj(const vtypeS& x, const vtypeS& y) const
{
	double ret = 1;
		
	for (size_t i = 0; i < D; ++i) {
		double tmp = 1.0;
		for (size_t j = 0; j < D; ++j)
			if (beta[i][j])
				tmp *= (1.0 - static_cast<double>(y[j]));
		
		ret *= (ck[i] * x[i] * tmp);  // prop fair utility
	}

	return std::log(ret);
}


template <size_t D>
bool
ALOHA<D>::constraints(const vtypeS& x, const vtypeS& y) const
{
	for (size_t i = 0; i < D; ++i)
	{
		double tmp = 1.0;

		for (size_t j = 0; j < D; ++j)
			if (beta[i][j])
				tmp *= (1.0 - static_cast<double>(x[j]));

		if (ck[i] * y[i] * tmp < Rmin[i])
			return false;
	}
	
	return true;
}
#endif
