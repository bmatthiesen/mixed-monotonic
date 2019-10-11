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


#ifndef _TIN_DM_H_
#define _TIN_DM_H_

#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>

#include "MMP.h"

template <size_t Dim>
class TIN_DM : public MMP<Dim>
{
	using typename MMP<Dim>::vtypeS;

	public:
		double alpha[Dim];
		std::array<double, Dim> beta[Dim];
		double sigma[Dim];

		TIN_DM() : MMP<Dim>() {}

	private:
		double MMPobj(const vtypeS& x, const vtypeS& y) const override;
};


template <size_t D>
double
TIN_DM<D>::MMPobj(const vtypeS& x, const vtypeS& y) const
{
	double ret = 1;

	for (size_t i = 0; i < D; ++i) {
		double tmp = alpha[i] * x[i] + std::inner_product(x.begin(), x.end(), beta[i].begin(), sigma[i]);
		tmp /= std::inner_product(y.begin(), y.end(), beta[i].begin(), sigma[i]);
		ret *= tmp;
	}

	return std::log2(ret);
}

#endif
