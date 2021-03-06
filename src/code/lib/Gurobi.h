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


#ifndef _GUROBI_H
#define _GUROBI_H

#include "gurobi_c++.h"

class GurobiEnv
{
	public:
		static size_t maxTries;
		static GRBEnv& getInstance();

		GurobiEnv() = delete;
		GurobiEnv(GurobiEnv const&) = delete;
		void operator=(GurobiEnv const&) = delete;

	private:
		static GRBEnv& _getInstance();
};

#endif
