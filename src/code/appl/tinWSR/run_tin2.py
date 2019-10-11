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
import tinPy
import tinDmPy
import tinPAPy
import tinRatespacePy
import tinMapelPy
import argparse
import os.path
from filelock import Lock

cf = h5py.File('../../../../data/chan.h5','r')

wp = int(os.getenv("SLURM_ARRAY_TASK_ID", "1003")) # 1003
savedir = os.getenv('JOB_HPC_SAVEDIR', "tmp")

if wp >= 900:
    Mapel = True
    wp = wp - 900
else:
    Mapel = False

idx = wp % 100
D = int(wp / 100) + 2

#parser = argparse.ArgumentParser()
#parser.add_argument("dimension", help="problem size", type=int)
#args = parser.parse_args()
#D = args.dimension

outfile = os.path.join(savedir, 'tin_dim{}_cidx{}_prec0.01.h5'.format(D,idx))
lockfile = outfile + '.lock'

beta = cf['channels'][idx,:D,:D]

dix = np.diag_indices(D)
alpha = beta[dix]
beta[dix] = 0

## MMP rel tol
#with h5py.File(outfile, 'r+') as f:
#    if 'mmp' in f and f['mmp']['Relative Tolerance'][0]:
#        f['mmp_reltol'] = f['mmp']
#        del f['mmp']
#
## MMP abs tol
#with h5py.File(outfile, 'r+') as f:
#    val = 'mmp' in f
#
#if not val:
#    w = getattr(tinPy, "TIN{}".format(D))(1e-2)
#    w.setPmax(1)
#    w.setChan(alpha, beta)
#    w.setRelTol(False)
#    w.optimize()
#
#    with h5py.File(outfile, 'r+') as f:
#        f.create_dataset('mmp', data = w.result())
#    del w
#
## DM rel tol
#with h5py.File(outfile, 'r+') as f:
#    if 'dm' in f and f['dm']['Relative Tolerance'][0]:
#        f['dm_reltol'] = f['dm']
#        del f['dm']
#
## DM abs tol
#with h5py.File(outfile, 'r+') as f:
#    val = 'dm' in f
#
#if not val:
#    g = getattr(tinDmPy, "TIN{}".format(D))(1e-2)
#    g.setPmax(1)
#    g.setChan(alpha, beta)
#    g.setRelTol(False)
#    g.optimize()
#
#    with h5py.File(outfile, 'r+') as f:
#        f.create_dataset('dm', data = g.result())
#    del g

if Mapel:
    # Mapel
    g = getattr(tinMapelPy, "TIN{}".format(D))(1e-2)
    g.setChan(1, alpha, beta)
    g.optimize()

    with Lock(lockfile):
        with h5py.File(outfile, 'r+') as f:
            if 'mapel' in f:
                del f['mapel']
            f.create_dataset('mapel', data = g.result())
    del g

else:
    # PA Ratespace
    g = getattr(tinRatespacePy, "TIN{}".format(D))(1e-2)
    g.setChan(1, alpha, beta)
    g.optimize()

    with Lock(lockfile):
        with h5py.File(outfile, 'r+') as f:
            if 'ratespace' in f:
                del f['ratespace']
            f.create_dataset('ratespace', data = g.result())
    del g
