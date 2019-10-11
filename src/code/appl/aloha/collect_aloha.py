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


import pathlib
import h5py

numchan = 100
p = pathlib.Path('/scratch/ws/bmatth-MMP/aloha')

with h5py.File(p / 'result.h5','w') as of:
    for fn in p.glob('*/**/*.h5'):
        try:
            with h5py.File(fn, 'r') as f:
                wpidx = f['input']['wpidx']
                for k in f.keys():
                    if k == 'input':
                        continue

                    if k not in of:
                        of.create_dataset(k, shape = (numchan, ), dtype = f[k].dtype)

                    of[k][wpidx] = f[k]
        except OSError:
            print('skipping %s' % fn)
