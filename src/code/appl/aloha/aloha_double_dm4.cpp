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

#include "aloha_double_dm.h"

int
main(int argc, char *argv[]) try
{
	const int Dim = 4;

	// parse input arguments
	if (argc != 4)
		throw std::invalid_argument("");

	const unsigned wpidx = std::stoi(argv[2]);
	const bool red = stob(argv[3]);
	const char *name = red ? "double_dm_red" : "double_dm";

	aloha_helper<ALOHA,Dim,2*Dim>(argv[1], wpidx, name, red);

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
	std::cerr << "Usage: " << argv[0] << " wpfile wpidx reduction" << std::endl;

	return -2;
}

