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
#include <thread>
#include <random>
#include <chrono>

#include "Gurobi.h"
#include "gurobi_c++.h"

GRBEnv&
GurobiEnv::getInstance()
{
	size_t tries = 0;

	std::mt19937_64 eng{std::random_device{}()};
    std::uniform_int_distribution<> dist{500, 10000}; // milliseconds

	while (true)
		try
		{
			return _getInstance();
		}
		catch (GRBException &e)
		{
			std::string msg = e.getMessage();

			if (e.getErrorCode() == GRB_ERROR_NO_LICENSE && msg.find("Failed to connect to token server") != std::string::npos && tries++ < maxTries)
			{
				auto sleep = dist(eng);

				std::cerr << "Gurobi Exception: " << msg << std::endl;
				std::cerr << "Retrying in: " << sleep * 1e-3 << "s" << "(try " << tries << " / " << maxTries << ")" << std::endl << std::endl;

				std::this_thread::sleep_for(std::chrono::milliseconds{sleep});
				continue;
			}

			throw;
		}
}

GRBEnv&
GurobiEnv::_getInstance()
{
	static GRBEnv instance;
	return instance;
}

size_t GurobiEnv::maxTries = 20;
