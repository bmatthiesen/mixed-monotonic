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

# cython: c_string_type=str, c_string_encoding=ascii
# vim: syntax=pyrex

from libcpp.vector cimport vector
from libcpp.string cimport string
from cython import boundscheck, wraparound
cimport libcpp
import numpy as np
cimport numpy as cnp
import h5py

cdef extern from "tin_ratespace.h":
    # https://stackoverflow.com/a/38567389/8372023
    #cdef enum Status "WSEE<DIM>::Status":
    #    Optimal, Unsolved, Infeasible

    # https://stackoverflow.com/a/41200186/8372023 
    cdef cppclass vtype`'`'DIM "std::array<double, DIM>":
        double& operator[](size_t)
        double* data()
        size_t size()

    cdef cppclass userdata`'`'DIM "userdata<DIM>":
        double N
        vtype`'`'DIM Pmax
        vtype`'`'DIM sigma
        vtype`'`'DIM alpha
        vtype`'`'DIM beta[DIM]
        void preprocess()

    cdef cppclass pa`'`'DIM "Tin_Ratespace<DIM>":
        pa`'`'DIM`'`'(void *)
        double getEpsilon()
        void setPrecision(double)
        void setUB(size_t, double)
        void optimize()

        vtype`'`'DIM xopt
        double optval
        unsigned long long iter
        unsigned long long lastUpdate
        const char* statusStr
        double runtime
        size_t max_queue_size
        const size_t data_size


cdef class TIN`'`'DIM:
    cdef pa`'`'DIM *_thisptr
    cdef userdata`'`'DIM *_ud
    cdef cnp.dtype result_dt

    @wraparound(False)
    @boundscheck(False)
    def __cinit__(self, double sigma):
        self._ud = new userdata`'`'DIM`'`'()
        self._thisptr = new pa`'`'DIM`'`'(self._ud)

        for i in range(DIM):
            self._ud.sigma[i] = sigma

        self.result_dt = np.dtype([
            ("Objective Value", np.float32),
            ("xopt", np.float32, DIM),
            ("Status", h5py.special_dtype(vlen=str)),
            ('Iterations', np.uint64),
            ('Last Update', np.uint64),
            ('Epsilon', np.float32),
            ('Runtime', np.float32),
            ('max_queue_size', np.uint64),
            ('data_size', np.uint16)])

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

        if self._ud != NULL:
            del self._ud

    @wraparound(False)
    @boundscheck(False)
    def setChan(self, P, double[:] alpha, double[:,::1] beta):
        if not (alpha.shape[0] == DIM and alpha.ndim == 1):
            raise ValueError('alpha must be 1D array with DIM elements')

        if not (beta.shape[0] == DIM and beta.shape[1] == DIM and beta.ndim == 2):
            raise ValueError("beta must be DIM()x`'DIM array")

        
        for i in range(DIM):
            self._ud.Pmax[i] = P
            self._ud.alpha[i] = alpha[i]

            for j in range(DIM):
                if i == j:
                    self._ud.beta[i][j] = 0.0
                else:
                    self._ud.beta[i][j] = beta[i][j]

            tmp = np.log2(1.0 + alpha[i] * P / self._ud.sigma[i])
            self._thisptr.setUB(i, tmp)
    
    def setPrecision(self, epsilon):
        self._thisptr.setPrecision(epsilon)

    def optimize(self):
        self._ud.preprocess();
        self._thisptr.optimize()

    def result(self):
        t = self._thisptr
        return np.array([(t.optval, self.xopt(), t.statusStr, t.iter, t.lastUpdate, t.getEpsilon(), t.runtime, t.max_queue_size, t.data_size)], dtype = self.result_dt)

    def getResultDt(self):
        return self.result_dt

    @wraparound(False)
    @boundscheck(False)
    def xopt(self):
        cdef size_t size = self._thisptr.xopt.size()
        cdef cnp.ndarray xopt = np.empty((size,))

        for i in range(size):
            xopt[i] = self._thisptr.xopt[i]

        return xopt
