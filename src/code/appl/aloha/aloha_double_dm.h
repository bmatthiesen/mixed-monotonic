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


#ifndef _ALOHA_DOUBLE_DM_H_
#define _ALOHA_DOUBLE_DM_H_

#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>

#include <iostream>
#include <sstream>

#include "MMP.h"
#include "aloha_helper.h"

using std::array;

template <size_t Dim>
class ALOHA : public MMPconstraints<2*Dim>
{
	using typename MMPconstraints<2*Dim>::vtypeS;
	using typename MMPconstraints<2*Dim>::RBox;

	public:
		double ck[Dim];
		double beta[Dim][Dim];
		double Rmin[Dim];
		
		ALOHA();

	private:
		double MMPobj(const vtypeS& x, const vtypeS& y) const override;
		bool constraints(const vtypeS& x, const vtypeS& y) const override;
		vtypeS feasiblePoint(const RBox& r) const override
			{ return r.lb(); }
};


template <size_t Dim>
ALOHA<Dim>::ALOHA()
	: MMPconstraints<2*Dim>()
{
	this->setUB(1);
	this->setLB(0);
}

template <size_t Dim>
double
ALOHA<Dim>::MMPobj(const vtypeS& x, const vtypeS&) const
{
	double ret = 1; // prop fair utility
		
	for (size_t i = 0; i < Dim; ++i) {
		double tmp = 1.0;
		for (size_t j = 0; j < Dim; ++j)
			if (beta[i][j] > 0)
				tmp *= x[Dim+j];
		
		ret *= (ck[i] * x[i] * tmp);  // prop fair utility
	}

	return std::log(ret);
}


template <size_t Dim>
bool
ALOHA<Dim>::constraints(const vtypeS& x, const vtypeS& y) const
{	
	for (size_t i = 0; i < Dim; ++i)
	{
		if (static_cast<double>(x[i]) + x[Dim+i] > 1)
			return false;
	}
	
	for (size_t i = 0; i < Dim; ++i)
	{
		double tmp = 1.0;
	
		for (size_t j = 0; j < Dim; ++j)
			if (beta[i][j] > 0)
				tmp *= y[Dim+j];
		
		if (ck[i] * y[i] * tmp < Rmin[i])
			return false;
	}

	return true;
}

#endif
