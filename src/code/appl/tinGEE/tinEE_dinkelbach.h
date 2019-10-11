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


#ifndef _TINEE_DINKELBACH_H_
#define _TINEE_DINKELBACH_H_

#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>

#include <iostream>

#include "MMP.h"

template <size_t Dim>
class TIN_Dinkelbach : public MMP<Dim>
{
	public:
		using typename MMP<Dim>::vtypeS;

		double alpha[Dim];
		std::array<double, Dim> beta[Dim];
		double sigma[Dim];

		double mu;
		const double outerTol = 1e-3;
		double Pc;

		unsigned long long outerIter;
		unsigned long long innerIter;

		TIN_Dinkelbach() : MMP<Dim>() {}

		//const typename H5::Result solve(const size_t idx);
		void solve(const size_t = 0);

	protected:
		double MMPobj(const vtypeS& x, const vtypeS& y) const override;

		using clk = std::chrono::high_resolution_clock;

		double lambda;
};


template <size_t D>
double
TIN_Dinkelbach<D>::MMPobj(const vtypeS& x, const vtypeS& y) const
{
	double ret = 1;

	for (size_t i = 0; i < D; ++i) {
		double tmp = alpha[i] * x[i] + std::inner_product(x.begin(), x.end(), beta[i].begin(), sigma[i]);
		tmp /= std::inner_product(y.begin(), y.end(), beta[i].begin(), sigma[i]);
		ret *= tmp;
	}

	ret = std::log2(ret);
	ret -= mu * lambda * std::accumulate(y.begin(), y.end(), static_cast<double>(0.0)); // the literal 0.0 is double, but let's be precise

	return ret;
}


template <size_t D>
void
TIN_Dinkelbach<D>::solve(const size_t)
{
	outerIter = 0;
	innerIter = 0;

	double cbv = 0;
	size_t queue_size = 0;

	std::chrono::time_point<clk> tic = clk::now();

	lambda = 0;

	do
	{
		++outerIter;
		this->optimize();

		cbv = this->optval - lambda * Pc;

		double numer = 1;
		for (size_t i = 0; i < D; ++i)
		{
			double tmp = alpha[i] * this->xopt[i];
			tmp /= std::inner_product(this->xopt.begin(), this->xopt.end(), beta[i].begin(), sigma[i]);
			tmp += 1;
			numer *= tmp;
		}
		numer = std::log2(numer);

		double denom = mu * std::accumulate(this->xopt.begin(), this->xopt.end(), 0.0) + Pc;

		lambda = numer / denom;

		innerIter += this->iter;
		queue_size = std::max(this->max_queue_size, queue_size);
	} while (std::abs(cbv) > outerTol);

	// update benchmark variables
	std::chrono::time_point<clk> toc = clk::now();
	this->max_queue_size = queue_size;
	this->runtime = std::chrono::duration<double>(toc - tic).count();

	/*H5::Result res(s.get());
	res.objective = lambda;
	res.iter = iter;
	res.lastUpdate = innerIter;
	res.runtime = std::chrono::duration<double>(toc - tic).count();*/

	// output
	std::cout << std::endl << std::endl;
	std::cout << "lambda: " << lambda << std::endl;
	std::cout << "iter: " << outerIter << std::endl;
	std::cout << "total inner iter: " << innerIter << std::endl;
	std::cout << "Time: " << this->runtime << std::endl;
}

#endif
