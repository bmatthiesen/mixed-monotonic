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
#include <stdexcept>

#include "aloha.h"
#include "aloha_helper.h"


int
main(int argc, char *argv[]) try
{
	const int Dim = 3;

	// parse input arguments
	if (argc != 3)
		throw std::invalid_argument("");

	const bool red = stob(argv[2]);
	const char *name = red ? "MMPred" : "MMP";

	for (unsigned wpidx = 0; wpidx < 100; ++wpidx)
		aloha_helper<ALOHA,Dim>(argv[1], wpidx, name, red);

	return 0;
}
catch (H5::Exception& e)
{
	std::cerr << "HDF Exception" << std::endl;
	e.printErrorStack();
	return -1;
}
catch(std::invalid_argument& e)
{
	std::cerr << "Usage: " << argv[0] << " wpfile reduction" << std::endl;

	return -2;
}
