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


## config

CXXFLAGS=-std=c++17
CPPFLAGS=-march=native -O3 -Wall -DNDEBUG -fpic -g -Wall #-Wno-unused-variable
#CPPFLAGS=-g  -Wall -DDEBUG -fpic
#CPPFLAGS=-ggdb -Wall -DDEBUG -fno-omit-frame-pointer
LDFLAGS=-g -march=native -m64
LDLIBS=-lstdc++ -lm

# Mosek
MSKFLAGS += -I${MSKHOME}/h -I$(abspath ${MSKHOME}/../../examples/c)
MSKLIBS += -L${MSKHOME}/bin -lmosek64

# Gurobi
GRBFLAGS += -I${GUROBI_HOME}/include

ifneq (,$(findstring taurus.hrsk.tu-dresden.de,$(shell hostname)))
  GRBLIBS = -L${HOME}/gurobi8/build
else
  GRBLIBS = -L${GUROBI_HOME}/src/build
endif

GRBLIBS +=-lgurobi_c++ -L${GUROBI_HOME}/lib -lgurobi81

# intel mkl
MKLLIBS +=  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl
MKLFLAGS +=  -m64 -I${MKLROOT}/include

.PHONY:
clean:
	rm -f *.o *.a
	rm -rf .d

### automatic dependency tracking
SRCS = $(wildcard *.cpp)

DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

COMPILE.c = $(CC) $(DEPFLAGS) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

%.o : %.c
%.o : %.c $(DEPDIR)/%.d
	$(COMPILE.c) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

%.o : %.cc
%.o : %.cc $(DEPDIR)/%.d
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

%.o : %.cxx
%.o : %.cxx $(DEPDIR)/%.d
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

%.o : %.cpp
%.o : %.cpp $(DEPDIR)/%.d
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d


include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS))))

