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
import tinFifoPy
import argparse
import os.path

cf = h5py.File('../../../../data/chan.h5','r')

wp = int(os.getenv("SLURM_ARRAY_TASK_ID", "5"))
savedir = os.getenv('JOB_HPC_SAVEDIR', "tmp")

idx = wp % 100
D = int(wp / 100) + 2

#parser = argparse.ArgumentParser()
#parser.add_argument("dimension", help="problem size", type=int)
#args = parser.parse_args()
#D = args.dimension

outfile = os.path.join(savedir, 'tin_dim{}_cidx{}_prec0.01.h5'.format(D,idx))

beta = cf['channels'][idx,:D,:D]

dix = np.diag_indices(D)
alpha = beta[dix]
beta[dix] = 0

with h5py.File(outfile, 'w') as f:
    f.create_dataset('dim', data = D)
    f.create_dataset('cidx', data = idx)

# MMP rel tol
w = getattr(tinPy, "TIN{}".format(D))(1e-2)
w.setPmax(1)
w.setChan(alpha, beta)
w.setRelTol(True)
w.optimize()

with h5py.File(outfile, 'r+') as f:
    f.create_dataset('mmp_reltol', data = w.result())
del w

# MMP abs tol
w = getattr(tinPy, "TIN{}".format(D))(1e-2)
w.setPmax(1)
w.setChan(alpha, beta)
w.setRelTol(False)
w.optimize()

with h5py.File(outfile, 'r+') as f:
    f.create_dataset('mmp', data = w.result())
del w

# MMP FiFo
w = getattr(tinFifoPy, "TIN{}".format(D))(1e-2)
w.setPmax(1)
w.setChan(alpha, beta)
w.setRelTol(False)
w.optimize()

with h5py.File(outfile, 'r+') as f:
    f.create_dataset('mmpFifo', data = w.result())
del w

# DM rel tol
g = getattr(tinDmPy, "TIN{}".format(D))(1e-2)
g.setPmax(1)
g.setChan(alpha, beta)
g.setRelTol(True)
g.optimize()

with h5py.File(outfile, 'r+') as f:
    f.create_dataset('dm_reltol', data = g.result())
del g

# DM abs tol
g = getattr(tinDmPy, "TIN{}".format(D))(1e-2)
g.setPmax(1)
g.setChan(alpha, beta)
g.setRelTol(False)
g.optimize()

with h5py.File(outfile, 'r+') as f:
    f.create_dataset('dm', data = g.result())
del g

# PA
p = getattr(tinPAPy, "TIN{}".format(D))(1e-2)
p.setPmax(1)
p.setChan(alpha, beta)
p.optimize()

with h5py.File(outfile, 'r+') as f:
    f.create_dataset('pa', data = p.result())
del p

# PA Ratespace
g = getattr(tinRatespacePy, "TIN{}".format(D))(1e-2)
g.setChan(1, alpha, beta)
g.optimize()

with h5py.File(outfile, 'r+') as f:
    f.create_dataset('ratespace', data = g.result())
del g

# Mapel
g = getattr(tinMapelPy, "TIN{}".format(D))(1e-2)
g.setChan(1, alpha, beta)
g.optimize()

with h5py.File(outfile, 'r+') as f:
    f.create_dataset('mapel', data = g.result())
del g
