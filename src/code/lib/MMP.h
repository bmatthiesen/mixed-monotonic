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


#ifndef _MMP_H
#define _MMP_H

#include <functional>
#include "BRB.h"

template <size_t Dim>
class _MMPbase : public BRB<Dim>
{
	public:
		using typename BRB<Dim>::vtype;
		using typename BRB<Dim>::RBox;
		using typename BRB<Dim>::PBox;
	    using BRB<Dim>::disableReduction;
		using vtypeS = typename PBox::vtype;

		_MMPbase() : BRB<Dim>() { disableReduction = true; }

	protected:
		virtual double MMPobj(const vtypeS& x, const vtypeS& y) const =0;

		// Return true, if feasible point is known. if feasible() == true, feasiblePoint() returns a feasible point.
		//bool feasible(const RBox& r) const override =0;
		//const vtype& feasiblePoint(const RBox& r) const override =0;

		// Return true if box is definitely empty.
		//bool isEmpty(const PBox& r) const override =0;

		void bound(RBox& r) const override final
			{ r.bound = MMPobj(r.ub(), r.lb()); }

		double obj(const RBox& r) const override final
			{ auto p = this->feasiblePoint(r); return MMPobj(p, p); }

		double red_alpha(const size_t, const double, const RBox&) const override
			{ throw std::runtime_error(ERR("not implemented")); } // TODO see commented out version below
		double red_beta(const size_t, const double, const RBox&) const override
			{ throw std::runtime_error(ERR("not implemented")); } // TODO see commented out version below
};

/* MMP objective, no constraints (except box) */
template <size_t Dim>
class MMP : public _MMPbase<Dim>
{
	public:
		using typename _MMPbase<Dim>::vtype;
		using typename _MMPbase<Dim>::RBox;
		using typename _MMPbase<Dim>::PBox;
		using typename _MMPbase<Dim>::vtypeS;

		using _MMPbase<Dim>::_MMPbase;

	protected:
		vtypeS feasiblePoint(const RBox& r) const override
			{	return r.lb(); }

		bool feasible(const RBox&) const override final
			{ return true; }

		bool isEmpty(const RBox&) const override final
			{ return false; }

		double red_alpha(const size_t i, const double gamma, const RBox& box) const override;
		double red_beta(const size_t i, const double gamma, const RBox& box) const override;
};

/* MMP objective, MMP constraints, no relaxation, infinite convergence */
template <size_t Dim>
class MMPconstraints : public _MMPbase<Dim>
{
	public:
		using typename _MMPbase<Dim>::vtype;
		using typename _MMPbase<Dim>::RBox;
		using typename _MMPbase<Dim>::PBox;
		using typename _MMPbase<Dim>::vtypeS;

		using _MMPbase<Dim>::_MMPbase;

	protected:
		// return true if G_i(x, y) <= 0 for all i
		//    where G_i are such that the constraints of the problem are "for all i and x: G_i(x, x) <= 0"
		virtual bool constraints(const vtypeS& x, const vtypeS& y) const =0;

		// this function should supply a candiate feasiblePoint. Its feasibility will be checked by constraints()
		//const vtype& feasiblePoint(const RBox& r) const override =0;


		bool feasible(const RBox& r) const override final
		{
			auto p = this->feasiblePoint(r);
			return constraints(p, p);
		}

		bool isEmpty(const RBox& r) const override final
		{
			return !constraints(r.lb(), r.ub());
		}

		double red_alpha(const size_t i, const double gamma, const RBox& box) const override;
		double red_beta(const size_t i, const double gamma, const RBox& box) const override;
};

template <class UnaryPredicate>
double
zero(UnaryPredicate feas)
{
	const double tol = 1e-1;
	double beta;
	double beta_min = 0;
	double beta_max = 1;
	bool ret;

	while (beta_max - beta_min > tol)
	{
		beta = (beta_max + beta_min) / 2;
		ret = feas(beta);

		if (ret)
			beta_min = beta;
		else
			beta_max = beta;
	}

	return beta_max;
}


template <size_t Dim>
double
MMP<Dim>::red_alpha(const size_t i, const double gamma, const RBox& box) const
{
	return zero([this, i, gamma, &box] (double alpha)
		{
			vtypeS p(box.ub());
			p[i] -= alpha * (static_cast<double>(box.ub(i)) - static_cast<double>(box.lb(i)));
			return this->MMPobj(p, box.lb()) >= gamma;
		});
}


template <size_t Dim>
double
MMP<Dim>::red_beta(const size_t i, const double gamma, const RBox& box) const
{
	return zero([this, i, gamma, &box] (double beta)
		{
			vtypeS p(box.lb());
			p[i] += beta * (static_cast<double>(box.ub(i))- static_cast<double>(box.lb(i)));
			return this->MMPobj(box.ub(), p) >= gamma;
		});
}

template <size_t Dim>
double
MMPconstraints<Dim>::red_alpha(const size_t i, const double gamma, const RBox& box) const
{
	return zero([this, i, gamma, &box] (double alpha)
		{
			vtypeS p(box.ub());
			p[i] -= alpha * (static_cast<double>(box.ub(i)) - static_cast<double>(box.lb(i)));
			return this->MMPobj(p, box.lb()) >= gamma && constraints(box.lb(), p);
		});
}

template <size_t Dim>
double
MMPconstraints<Dim>::red_beta(const size_t i, const double gamma, const RBox& box) const
{
	return zero([this, i, gamma, &box] (double beta)
		{
			vtypeS p(box.lb());
			p[i] += beta * (static_cast<double>(box.ub(i)) - static_cast<double>(box.lb(i)));
			return this->MMPobj(box.ub(), p) >= gamma && constraints(p, box.ub());
		});
}

#endif
