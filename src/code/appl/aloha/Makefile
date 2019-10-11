# Copyright (C) 2018-2019 Bho Matthiesen, Christoph Hellings
# 
# This program is used in the article:
#
# Bho Matthiesen, Christoph Hellings, Eduard A. Jorswieck, and Wolfgang
# Utschick, "Mixed Monotonic Programming for Fast Global Optimization,"
# submitted to IEEE  Transactions on Signal Processing.
# 
# 
# License:
# This program is licensed under the GPLv2 license. If you in any way use this
# code for research that results in publications, please cite our original
# article listed above.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


EXE := aloha aloha4 aloha_double_dm aloha_double_dm4 aloha_gp aloha_gp4 aloha_double
include ../common.mk

CXX=h5c++
CC=h5c++

dgopt.o scopt-ext.o: CC=cc

# Mosek
vpath scopt-ext.c $(abspath ${MSKHOME}/../../examples/c)
scopt-ext.o dgopt.o: CPPFLAGS:=$(filter-out -Wall -DDEBUG,$(CPPFLAGS))
scopt-ext.o dgopt.o: CPPFLAGS+=$(MSKFLAGS) -Wno-unused-result

aloha_gp aloha_gp4: LDLIBS+=$(MSKLIBS)
aloha_gp aloha_gp4: Mosek.o scopt-ext.o dgopt.o
aloha_gp.o aloha_gp4.o Mosek.o: CPPFLAGS+=$(MSKFLAGS)

aloha_double: LDLIBS+=$(MKLLIBS)
aloha_double.o: CPPFLAGS+=$(MKLFLAGS)
