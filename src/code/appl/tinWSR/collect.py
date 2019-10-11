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
import numpy as np

ofn = 'result.h5'
#p = pathlib.Path('/scratch/p_mwrc/diss/benchPA')
p = pathlib.Path('tmp')

nchan = 100
maxdim = 10

tin_algos = ('dm', 'dm_reltol', 'mmp', 'mmp_reltol', 'pa', 'ratespace', 'mapel', 'mmpFifo')
gee_algos = ('gee_mmp', 'gee_dinkelbach')

of = h5py.File(p / ofn, 'w')
for alg in tin_algos + gee_algos:
    tmp = of.create_group(alg)

    tmp.create_group('raw_results')

    rt = tmp.create_dataset('runtime', dtype='f4', shape = (maxdim, nchan), fillvalue = np.nan)
    mem = tmp.create_dataset('memory', dtype='f16', shape = (maxdim, nchan), fillvalue = np.nan)
    it = tmp.create_dataset('iterations', dtype='f8', shape = (maxdim, nchan), fillvalue = np.nan)


def collect(glob, algos):
    for fn in p.glob(glob):
        with h5py.File(fn,'r') as f:

            d = int(f['dim'][...])
            cidx = int(f['cidx'][...])

            for alg in algos:
                try:
                    if str(d) not in of[alg]['raw_results']:
                        of[alg]['raw_results'].create_dataset(str(d), shape = (nchan,), dtype = f[alg].dtype)

                    of['{}/raw_results/{}'.format(alg, d)][cidx] = f[alg][...]

                    of[alg]['runtime'][d-1, cidx] = f[alg]['Runtime'][...]
                    of[alg]['memory'][d-1, cidx] = f[alg]['max_queue_size'][...] * f[alg]['data_size'][...]
                    of[alg]['iterations'][d-1, cidx] = f[alg]['Iterations'][...]

                except KeyError:
                    pass

collect('tin_dim*h5', tin_algos)
collect('gee_dim*h5', gee_algos)

of.close()
