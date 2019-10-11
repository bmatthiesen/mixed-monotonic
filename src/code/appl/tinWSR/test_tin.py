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


import numpy as np
import tinPy
import tinFifoPy
import tinDmPy
import tinPAPy
import tinRatespacePy
import tinMapelPy

beta = np.array([[8.3401758e+02, 5.9968562e+00, 9.5184622e+00, 6.0737956e-01],
[1.3587096e+00, 3.9182301e+01, 2.0014184e-02, 1.6249435e+00],
[3.8521406e-01, 4.6761915e-01, 8.7457578e+03, 1.8704400e+00],
[1.2729254e-01, 2.1447293e-02, 3.1017335e-02, 1.2471862e+02]])

beta = np.array([[8.3401758e+02, 5.9968562e+00, 9.5184622e+00],
[1.3587096e+00, 3.9182301e+01, 2.0014184e-02],
[3.8521406e-01, 4.6761915e-01, 8.7457578e+03]
])

dix = np.diag_indices(3)
alpha = beta[dix]
beta[dix] = 0

#alpha = np.array([5.0091360e-0, 2.4832193e+00, 2.6642494e+00])
#beta = np.array([[0, 7.9828486e-01, 1.0228757e+00], [1.7795534e-01, 0, 8.9789916e-02], [6.0391173e-01, 7.3569861e-01, 0]])

w = tinPy.TIN3(1e-2)
w.setPmax(1)
w.setChan(alpha, beta)
w.setPrecision(1e-1)
w.optimize()

#g = tinRatespacePy.TIN3(1e-2)
#g.setChan(1, alpha, beta)
#g.setPrecision(1e-4)
#g.optimize()

#g = tinMapelPy.TIN3(1e-2)
#g.setChan(1, alpha, beta)
#g.setPrecision(1e-4)
#g.optimize()

#q = tinDmPy.TIN4(1e-2)
#q.setPmax(1)
#q.setChan(alpha, beta)
#q.optimize()

g = tinPAPy.TIN3(1e-2)
g.setPmax(1)
g.setChan(alpha, beta)
g.setPrecision(1e-1)
g.optimize()
