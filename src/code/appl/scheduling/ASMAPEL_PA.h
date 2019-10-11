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


#ifndef _ASMAPEL_PA_H
#define _ASMAPEL_PA_H

#include <iostream>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <array>
#include <utility>
#include <iterator>
#include <cstdio>
#include <chrono>

#include "util.h"


/* base types */

namespace smapel {
typedef double basetype;
enum PAstatus { OPTIMAL, UNSOLVED, MAXITER };

/* main class */
template <size_t Dim>
class PA
{
	public:
		using vtype = std::array<basetype, Dim>;

		PA(void *ud = nullptr)
			: output(true), output_every(100), MaxIter(-1), user_data(ud), xopt {}, data_size((Dim+1)*sizeof(basetype)), epsilon(1e-2), proj_delta(1e-10), tol(0.001), ub {}, obj(NULL), inG(NULL), inH(NULL), doshift(false), xshift {}
		{ setStatus(UNSOLVED); }

		virtual ~PA() {};

		struct Vertex
		{
			vtype v;
			basetype obj;

			Vertex() {};

			Vertex(const vtype &vert, basetype b = -std::numeric_limits<double>::infinity())
				: v(vert), obj(b)
			{}

			Vertex(vtype &&vert, basetype b = -std::numeric_limits<double>::infinity())
				: v(std::move(vert)), obj(b)
			{ }
		};

		using VertexSet = std::vector<Vertex>;

		using PAfun = basetype (*)(const vtype &, void *);
		using PAboolfun = bool (*)(const vtype &, void *);

		// vector length
		constexpr size_t dim() { return Dim; };

		// parameter setter
		void setPrecision(basetype epsilon)
			{ this->epsilon = epsilon; setStatus(UNSOLVED); }

		void setObjective(PAfun fp)
			{ obj = fp; setStatus(UNSOLVED); }

		void setInG(PAboolfun fp)
			{ inG = fp; setStatus(UNSOLVED); }

		void setInH(PAboolfun fp)
			{ inH = fp; setStatus(UNSOLVED); }

		void setUB(vtype v)
			{ub = v; setStatus(UNSOLVED); }

		void setUB(basetype e)
		{
			for (auto &b : ub)
				b = e;

			setStatus(UNSOLVED);
		}

		void setTol(double t)
		{
			tol = t;
			setStatus(UNSOLVED);
		}

		void setUB(size_t idx, basetype e)
			{ub[idx] = e; setStatus(UNSOLVED); }


		void setUserData(void *ud)
			{ user_data = ud; setStatus(UNSOLVED); }

		void setShift(basetype e);

		void setShift(vtype v = {})
		{ xshift = v; doshift = true; }

		void unsetShift()
			{ doshift = false; }

		// parameter
		bool output;
		unsigned long long output_every;
		unsigned long long MaxIter;

		void *user_data;

		// result
		vtype xopt;
		basetype optval;
		unsigned long long iter;
		unsigned long long lastUpdate;
		PAstatus status;
		const char *statusStr;
		double runtime; // in seconds
		size_t max_queue_size;
		const size_t data_size;

		// run PA algorithm
		void optimize(vtype startPoint)
			{ optimize(true, std::forward(startPoint)); }

		void optimize()
			{ optimize(false, {}); }

		void printResult();
		double getEpsilon() const { return epsilon; }

	protected:
		// parameter
		basetype epsilon;
		const double proj_delta;
		double tol;

		vtype ub;
		PAfun obj;
		PAboolfun inG, inH;

		bool doshift;
		vtype xshift;

		typedef std::chrono::high_resolution_clock clock;

		// functions
		struct projection_return
		{
			vtype ret;
			vtype cand;
		};
		virtual projection_return projection(const vtype& x, void*) =0;
		void optimize(bool startPointPresent, typename PA<Dim>::vtype startPoint);

		vtype shift(const vtype &x) const
		{
			if (!doshift)
				return x;

			vtype res;
			std::transform(
				x.begin(),
				x.end(),
				xshift.begin(),
				res.begin(),
				[this] (const basetype &a, const basetype &b)
				{
					return a - b;
				}
			);

			return res;
		}

		void setStatus(const PAstatus s)
		{
			status = s;

			switch (s)
			{
				case OPTIMAL:
					statusStr = "Optimal";
					break;

				case UNSOLVED:
					statusStr = "Undefined";
					break;

				case MAXITER:
					statusStr = "Max Iterations";
			}
		}
};

#include "bits/ASMAPEL_PA.cpp"
}
#endif
