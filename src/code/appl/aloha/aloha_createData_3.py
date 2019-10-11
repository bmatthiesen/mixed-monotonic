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


import h5py
import numpy as np
import itertools as it

import sys
sys.path.append('python')
import alohaPy

dim = 3;
numchan = 100;
interference_level = 0
SNR = 10
sigma = .05

def crandn(*args, **kwargs):
    return 1/np.sqrt(2) * (np.random.randn(*args, **kwargs) + 1j * np.random.randn(*args, **kwargs))

with h5py.File('../../data/aloha' + str(dim) + '.h5', 'w') as ds:
    ck = np.log2(1 + np.abs(crandn(numchan, dim))**2 * SNR)
    ds.create_dataset("ck", data = ck, maxshape = (None, dim))

    # TODO keine leeren zeilen zulassen
    beta = np.random.rand(numchan,dim,dim) > interference_level
    idx = np.diag_indices(dim)
    for i in range(beta.shape[0]):
        beta[i][idx] = False
    ds.create_dataset("beta", data = beta, maxshape = (None, dim, dim))

    # choose Rmin around 4/27 of ck with standard deviation .05
    Rmin = (4/27 + np.random.randn(numchan, dim) * sigma) * ck

    btmp = np.asarray(beta, dtype=np.dtype('i'))
    a = alohaPy.ALOHA3()
    for i in range(numchan):
        while True:
            a.setParam(ck[i], btmp[i], Rmin[i])
            
            if a.checkFeas():
                break

            print("updating %d" % i)
            Rmin[i] = (4/27 + np.random.randn(dim) * sigma) * ck[i]

    ds.create_dataset("Rmin", data = Rmin, maxshape = (None, dim))
