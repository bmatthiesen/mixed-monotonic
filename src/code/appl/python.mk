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


PY=/usr/bin/python3
PYEXT=cpython-37m-x86_64-linux-gnu.so

EXTNAMES=$(addsuffix $(PYEXT),$(EXT))

.PHONY: python
python: $(addsuffix .pyx,$(EXT)) setup.py
	$(PY) setup.py build_ext --inplace

%Py.pyx: %Py.pyx.m4 %Py.2.pyx.in %Py.3.pyx.in %Py.4.pyx.in %Py.5.pyx.in %Py.6.pyx.in %Py.7.pyx.in %Py.8.pyx.in %Py.9.pyx.in %Py.10.pyx.in %Py.11.pyx.in %Py.12.pyx.in %Py.13.pyx.in %Py.14.pyx.in %Py.15.pyx.in %Py.16.pyx.in %Py.17.pyx.in %Py.18.pyx.in %Py.19.pyx.in %Py.20.pyx.in
	cat $(wordlist 2,$(words $^),$^) > $@
	@rm $(wordlist 2,$(words $^),$^)

%.pyx.in: $(addsuffix .pyx.m4,$(basename $*))
	m4 -P --define=DIM=$(subst .,,$(suffix $*)) $(addsuffix .pyx.m4,$(basename $*)) > $@

pyclean:
	rm -f *.pyx *.$(PYEXT) $(addsuffix .cpp,$(EXT))

