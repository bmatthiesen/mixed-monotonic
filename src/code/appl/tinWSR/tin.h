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


#ifndef _TIN_H_
#define _TIN_H_

#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>

#include <iostream>

#include "MMP.h"

template <size_t Dim, bool propFair = false>
class TIN : public MMP<Dim>
{
	using typename MMP<Dim>::vtypeS;

	public:
		double alpha[Dim];
		std::array<double, Dim> beta[Dim];
		double sigma[Dim];

		TIN() : MMP<Dim>() {}

	private:
		double MMPobj(const vtypeS& x, const vtypeS& y) const override;
};


template <size_t D, bool propFair>
double
TIN<D,propFair>::MMPobj(const vtypeS& x, const vtypeS& y) const
{
	double ret = 1;

	for (size_t i = 0; i < D; ++i) {
		double tmp = 1 + alpha[i] * x[i] / (std::inner_product(y.begin(), y.end(), beta[i].begin(), sigma[i]) + beta[i][i]*(static_cast<double>(x[i])-y[i]));

		if constexpr(propFair)
			ret *= std::log2(tmp);
		else
			ret *= tmp;
	}

	if constexpr(propFair)
		return ret;
	else
		return std::log2(ret);
}

#endif
