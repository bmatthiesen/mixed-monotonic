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


#ifndef _BRB_H
#define _BRB_H

#include <vector>
#include <array>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cstdio>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <clocale>
#include <queue>
#include <memory>

#include "util.h"

using std::cout;
using std::endl;

#define ERR(s)  std::string(__func__) + ": " + (s)

template <typename T>
std::ostream& operator<< (std::ostream &out, const std::vector<T>& v)
{
	out << "[";
	for (auto& e : v)
		out << " " << e;
	out << " ]";

	return out;
}


template<typename BRB>
struct compare_RBox
{
    bool operator()(const typename BRB::RBox& n1, const typename BRB::RBox& n2) const
    { return n1.bound < n2.bound; }
};

template <size_t Dim, typename storageT = float>
class BRB
{
	public:
    using vtype = std::array<double, Dim>;

    struct PBox
    {
		using vtype = std::array<storageT, Dim>;
		vtype lb, ub;
    };

    class RBox
    {
        std::unique_ptr<PBox> data_;
    public:
        RBox() : data_(std::make_unique<PBox>()) {};
        RBox(std::unique_ptr<PBox>&& data) : data_(std::move(data)) {};

        std::unique_ptr<PBox> move_data() {
			return std::move(data_);
		}

        storageT lb(size_t index) const { return data_->lb[index]; }
        storageT& lb(size_t index) { return data_->lb[index]; }
        const typename PBox::vtype& lb() const { return data_->lb; }
		typename PBox::vtype& lb() { return data_->lb; }

        storageT ub(size_t index) const { return data_->ub[index]; }
        storageT& ub(size_t index) { return data_->ub[index]; }
        const typename PBox::vtype& ub() const { return data_->ub; }
		typename PBox::vtype& ub() { return data_->ub; }

        storageT bound;
    };

#ifdef FIFO
	using RType = std::queue<RBox, std::deque<RBox>>;
#else
    using RType = std::priority_queue<RBox, std::vector<RBox>, compare_RBox<BRB>>;
#endif


    class PType
    {
    public:
        PType() : len(1) { };

        PBox& operator[](const size_t index)
        {
            if (len == 1)
                return P1;
            else
                return P2[index];
        }

        const PBox& operator[](const size_t index) const
        {
            if (len == 1)
                return P1;
            else
                return P2[index];
        }

        void use1() { len = 1; }
        void use2() { len = 2; }
        size_t size() { return len; }

    private:
        PBox P1;
        std::array<PBox, 2> P2;
        size_t len;
    };

    public:
		BRB();

		virtual ~BRB() {};

		// types
		enum class Status { Optimal, Unsolved, Infeasible };

		// parameter setter
		void setPrecision(const double eta);
		void setLB(const vtype& v);
		void setLB(const double e);
		void setLB(const size_t idx, const double e);
		void setUB(const vtype& v);
		void setUB(const double e);
		void setUB(const size_t idx, const double e);

		// parameter
		bool output;
		unsigned long long outputEvery;
		bool disableReduction;
		bool useRelTol;
		bool enablePruning;

		// result
		vtype xopt;
		double optval;
		unsigned long long iter, lastUpdate;
		Status status;
		const char *statusStr;
		double runtime; // in seconds
		size_t max_queue_size;
		const size_t data_size;

		// run algorithm
		virtual void optimize(bool startFromXopt = false);
		virtual void printResult() const;

		constexpr size_t dim() const { return Dim; }
		double getEpsilon() const { return epsilon; }

	protected:
		// parameter
		typedef std::chrono::high_resolution_clock clock;
		std::chrono::time_point<clock> tic;
		vtype lb, ub;
		double epsilon;

		// functions
		void setStatus(const Status s);

		void reduction(RBox& red, const double gamma) const; // false if reduced box is empty

		virtual void bound(RBox& r) const =0;

		virtual double red_alpha(const size_t i, const double gamma, const RBox& box) const =0;
		virtual double red_beta(const size_t i, const double gamma, const RBox& box) const =0;

		virtual typename PBox::vtype feasiblePoint(const RBox& r) const =0;

		virtual bool isEmpty(const RBox& r) const =0;
		virtual bool feasible(const RBox& r) const =0;
		virtual double obj(const RBox& r) const =0;

		virtual void checkpoint() const;

		void prune(RType& R);
};


template <size_t Dim, typename storageT>
BRB<Dim, storageT>::BRB() : output(true), outputEvery(1e6), disableReduction(false), useRelTol(true), enablePruning(false), data_size((2*Dim+1) * sizeof(storageT)), epsilon(1e-2)
{
	setStatus(Status::Unsolved);
	std::setlocale(LC_NUMERIC, "en_US.UTF-8");
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::setPrecision(const double eps)
{
	epsilon = eps;
	setStatus(Status::Unsolved);
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::setLB(const vtype& v)
{
	lb = v;
	setStatus(Status::Unsolved);
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::setLB(const double e)
{
	for (auto &b : lb)
		b = e;

	setStatus(Status::Unsolved);
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::setLB(const size_t idx, const double e)
{
	lb[idx] = e;
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::setUB(const vtype& v)
{
	ub = v;
	setStatus(Status::Unsolved);
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::setUB(const double e)
{
	for (auto &b : ub)
		b = e;

	setStatus(Status::Unsolved);
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::setUB(const size_t idx, const double e)
{
	ub[idx] = e;
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::printResult() const
{
	cout << "Status: " << statusStr << endl;
	cout << "Optval: " << optval << endl;

	cout << "X*: [";
	std::for_each(xopt.begin(), xopt.end(), [] (const double &a) { cout << " " << a; });
	cout << " ]" << endl;

	cout << "Iter: " << iter << endl;
	cout << "Solution found in iter: " << lastUpdate << endl;

	cout << "Runtime: " << runtime << " sec" << endl;
}

inline
double
calcTolerance(bool u, double o, double e)
{
	return u ? (o != 0 ? (1+e)*o : o+std::numeric_limits<double>::epsilon()) : o+e;
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::prune(RType& R)
{
	RType Rnew;
	const double gamma = calcTolerance(useRelTol, optval, epsilon);

	// prune if bound < gamma

#ifdef FIFO
	while (!R.empty() && R.front().bound >= gamma)
	{
		Rnew.push(std::move(const_cast<RBox&>(R.front())));
#else
	while (!R.empty() && R.top().bound >= gamma)
	{
	    // Supposedly this is legal: https://stackoverflow.com/a/20149745/620382
		Rnew.push(std::move(const_cast<RBox&>(R.top())));
#endif
		R.pop();
	}

	std::swap(R, Rnew);

	return;
}

/*
 * Algorithm BRB in Section 7.6 of Tuy et al., "Monotonic Optimization: Branch
 * and Cut Methods", in: Audet et al., "Essays and Surveys in Global
 * Optimization," Springer 2005.
 */
template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::optimize(bool startFromXopt)
{
	PType P;
	RType R;
	bool doPrune = false;

	if (startFromXopt)
	{
		throw std::runtime_error(ERR("not implemented"));

		// set from CP: runtime, xopt, iter
		//tic = clock::now() - std::chrono::duration_cast<clock::duration>(std::chrono::duration<double>(runtime));
		//optval = obj(xopt);
		//lastUpdate = iter;

	}
	else
	{
		tic = clock::now();
		iter = lastUpdate = 0;
		optval = -std::numeric_limits<double>::infinity();
	}

	setStatus(Status::Unsolved);

	// step 0
	P.use1();
	std::copy(lb.begin(), lb.end(), P[0].lb.begin());
	std::copy(ub.begin(), ub.end(), P[0].ub.begin());

	MiniPool<PBox> pool;
    while (true)
	{
		iter++;

		// step 1: reduce & bound
		for (size_t i = 0; i < P.size(); i++)
		{
			RBox box(pool.get());

			box.lb() = P[i].lb;
			box.ub() = P[i].ub;

			bound(box); // set bound
			const double gamma = calcTolerance(useRelTol, optval, epsilon);

			if (box.bound < gamma || isEmpty(box)) {
				pool.put(box.move_data());
				continue; // skip boxes containing no feasible points
			}

			if (!disableReduction)
			{
				reduction(box, gamma); // update lb, ub
				bound(box); // TODO calculating bound twice!
			}

			// step 2: update CBV, update R
			if (feasible(box))
			{
				double tmp = obj(box);

				if (tmp > optval)
				{
					optval = tmp;
					auto tmp = feasiblePoint(box);
					std::copy(tmp.begin(), tmp.end(), xopt.begin());

					doPrune = enablePruning;
					lastUpdate = iter;
				}
			}

			R.push(std::move(box));
		}

		if (doPrune && iter - lastUpdate > 10000)
		{
			doPrune = false;

			std::printf("PRUNE: %'zu", R.size());
			prune(R);
			std::printf(" %'zu\n", R.size());
		}

		// step 3: terminate
#ifdef FIFO
		if (R.empty())
#else
		if (R.empty() || R.top().bound - calcTolerance(useRelTol, optval, epsilon) < 0.0) // TODO calcTolerance() correct?
#endif
		{
			if (optval == -std::numeric_limits<double>::infinity())
				setStatus(Status::Infeasible);
			else
				setStatus(Status::Optimal);

			break;
		}

		// step 4: select box & branch
		{
#ifdef FIFO
			typename RType::const_reference M = R.front(); // argmax
#else
			typename RType::const_reference M = R.top(); // argmax
#endif

			if (output && iter % outputEvery == 0)
				std::printf("%'8llu  %'8lu  %11g  %11g  (%11g | %'8llu) | Peak RSS: %'zu\n", iter, R.size(), optval, M.bound, (M.bound - optval), lastUpdate, getPeakRSS());


			vtype tmp;
			std::transform(M.ub().begin(), M.ub().end(), M.lb().begin(), tmp.begin(), [] (double x, double y) { return x-y; });
			auto max = std::max_element(tmp.begin(), tmp.end());
			size_t jk = std::distance(tmp.begin(), max);

			auto vk = M.lb(jk) + static_cast<double>(*max) / 2;

			P.use2();

			P[0].lb = M.lb();
			P[0].ub = M.ub();
			P[0].ub[jk] = vk;

			P[1].lb = M.lb();
			P[1].lb[jk] = vk;
			P[1].ub = M.ub();

			pool.put(const_cast<RBox&>(M).move_data());
			R.pop();
		}
	}

	runtime = std::chrono::duration<double>(clock::now() - tic).count();
	max_queue_size = pool.size() + R.size();

	if (output)
	{
		std::printf("%'8llu  %'8lu  %11g  %11g  (%11g | %'8llu) | Peak RSS: %'zu\n", iter, R.size(), optval, std::nan("1"), 0.0, lastUpdate, getPeakRSS());
		std::cout << "max_queue_size = " << max_queue_size << std::endl;
		printResult();
	}
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::reduction(RBox& red, const double gamma) const
{
	// compute red.lb
	for (size_t i = 0; i < Dim; ++i)
	{
		const double alpha = red_alpha(i, gamma, red);
		red.lb(i) = red.ub(i) - alpha * (red.ub(i) - red.lb(i));
	}

	// compute red.ub
	for (size_t i = 0; i < Dim; ++i)
	{
		const double beta = red_beta(i, gamma, red);
		red.ub(i) = red.lb(i) + beta * (red.ub(i) - red.lb(i));
	}
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::setStatus(const Status s)
{
	status = s;

	switch (s)
	{
		case Status::Optimal:
			statusStr = "Optimal";
			break;

		case Status::Unsolved:
			statusStr = "Unsolved";
			break;

		case Status::Infeasible:
			statusStr = "Infeasible";
			break;
	}
}

template <size_t Dim, typename storageT>
void
BRB<Dim, storageT>::checkpoint() const
{
	return;
}


#endif
