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
import geePy
import dinkelbachPy
import argparse
import os.path

cf = h5py.File('../../../../data/chan.h5','r')

wp = int(os.getenv("SLURM_ARRAY_TASK_ID", "206"))
savedir = os.getenv('JOB_HPC_SAVEDIR', "tmp")

idx = wp % 100
D = int(wp / 100) + 2

outfile = os.path.join(savedir, 'gee_dim{}_cidx{}_prec0.01.h5'.format(D,idx))

beta = cf['channels'][idx,:D,:D]

dix = np.diag_indices(D)
alpha = beta[dix]
beta[dix] = 0

with h5py.File(outfile, 'w') as f:
    f.create_dataset('dim', data = D)
    f.create_dataset('cidx', data = idx)

# GEE MMP abs tol
w = getattr(geePy, "GEE_TIN{}".format(D))(5,1,1e-3)
w.setPmax(1)
w.setChan(alpha, beta)
w.setRelTol(False)
w.optimize()

with h5py.File(outfile, 'r+') as f:
    f.create_dataset('gee_mmp', data = w.result())
del w

# GEE Dinkelbach abs tol
w = getattr(dinkelbachPy, "Tin_Dinkelbach{}".format(D))(15,1,1)
w.setPmax(1)
w.setChan(alpha, beta)
w.setRelTol(False)
w.optimize()

with h5py.File(outfile, 'r+') as f:
    f.create_dataset('gee_dinkelbach', data = w.result())
del w
