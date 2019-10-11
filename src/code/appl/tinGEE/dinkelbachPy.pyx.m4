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

cdef extern from "tinEE_dinkelbach.h":
    # https://stackoverflow.com/a/38567389/8372023
    #cdef enum Status "WSEE<DIM>::Status":
    #    Optimal, Unsolved, Infeasible

    # https://stackoverflow.com/a/41200186/8372023 
    cdef cppclass vtype "std::array<double, DIM>":
        double& operator[](size_t)
        double* data()
        size_t size()

    cdef cppclass tinEE_dm`'`'DIM "TIN_Dinkelbach<DIM>":
        double mu
        double Pc
        double alpha[DIM]
        vtype beta[DIM]
        double sigma[DIM]

        libcpp.bool useRelTol
        libcpp.bool disableReduction

        vtype xopt
        double optval
        unsigned long long innerIter
        unsigned long long outerIter
        unsigned long long lastUpdate
        const char* statusStr
        double runtime
        size_t max_queue_size
        const size_t data_size

        gee_tin`'`'DIM`'`'()
        void setUB(double)
        void setLB(double)
        void setPrecision(double)
        void solve()
        double getEpsilon()


cdef class Tin_Dinkelbach`'`'DIM:
    cdef tinEE_dm`'`'DIM *_thisptr
    cdef cnp.dtype result_dt

    @wraparound(False)
    @boundscheck(False)
    def __cinit__(self, double mu, double psi, double sigma):
        self._thisptr = new tinEE_dm`'`'DIM`'`'()

        self._thisptr.setLB(0)
        self._thisptr.setPrecision(1e-2)
        self._thisptr.useRelTol = True
        self._thisptr.disableReduction = True

        self._thisptr.Pc = psi
        self._thisptr.mu = mu

        for i in range(DIM):
            self._thisptr.sigma[i] = sigma

        self.result_dt = np.dtype([
            ("Objective Value", np.float32),
            ("xopt", np.float32, DIM),
            ("Status", h5py.special_dtype(vlen=str)),
            ('Iterations', np.uint64),
            ('Dinkelbach Iterations', np.uint64),
            ('Epsilon', np.float32),
            ('Relative Tolerance', np.bool),
            ('Runtime', np.float32),
            ('max_queue_size', np.uint64),
            ('data_size', np.uint16)])

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    @wraparound(False)
    @boundscheck(False)
    def setChan(self, double[:] alpha, double[:,::1] beta):
        if not (alpha.shape[0] == DIM and alpha.ndim == 1):
            raise ValueError('alpha must be 1D array with DIM elements')

        if not (beta.shape[0] == DIM and beta.shape[1] == DIM and beta.ndim == 2):
            raise ValueError("beta must be DIM()x`'DIM array")

        
        for i in range(DIM):
            self._thisptr.alpha[i] = alpha[i]

            for j in range(DIM):
                self._thisptr.beta[i][j] = beta[i][j]

    def setPmax(self, P):
        self._thisptr.setUB(P)
    
    def setPrecision(self, epsilon):
        self._thisptr.setPrecision(epsilon)

    def setRelTol(self, val):
        self._thisptr.useRelTol = val

    def optimize(self):
        self._thisptr.solve()

    def result(self):
        t = self._thisptr
        return np.array([(t.optval, self.xopt(), t.statusStr, t.innerIter, t.outerIter, t.getEpsilon(), t.useRelTol, t.runtime, t.max_queue_size, t.data_size)], dtype = self.result_dt)

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
