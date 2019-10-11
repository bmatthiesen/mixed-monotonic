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


#include <iostream>
#include "Mosek.h"

void MSKAPI
printstr(void *, const char str[])
{
	std::cout << str;
}

void
mskcall(const MSKrescodee res)
{
	char buf[MSK_MAX_STR_LEN];

	if (res != MSK_RES_OK)
	{
		MSK_getcodedesc(res, buf, NULL);

		throw MSK_error(res, buf);
	}
}


MSK_error::MSK_error(const MSKrescodee err, const std::string str)
	: std::runtime_error("Mosek error: " + str + " (" + std::to_string(err) + ")"), rescode(err)
{}


MosekEnv::MosekEnv() try
	: env(nullptr)
{
	mskcall(MSK_makeenv(&env, nullptr));
}
catch(...)
{
	destroy();
}

MosekEnv::~MosekEnv() noexcept
{
	destroy();
}

void
MosekEnv::destroy() noexcept
{
	if (env != nullptr)
		MSK_deleteenv(&env);
}

MSKenv_t
MosekEnv::getEnv()
{
	static MosekEnv e;

	return e.env;
}
