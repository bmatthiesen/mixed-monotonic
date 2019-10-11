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
m4_define(DIM2, m4_incr(DIM))

from libcpp.vector cimport vector
from libcpp.string cimport string
from cython import boundscheck, wraparound
cimport libcpp
import numpy as np
cimport numpy as cnp
import h5py

cdef extern from "tin_PA.h":
    # https://stackoverflow.com/a/38567389/8372023
    #cdef enum Status "WSEE<DIM>::Status":
    #    Optimal, Unsolved, Infeasible

    # https://stackoverflow.com/a/41200186/8372023 
    cdef cppclass vtype`'`'DIM "std::array<double, DIM>":
        double& operator[](size_t)
        double* data()
        size_t size()
    cdef cppclass vtype`'`'DIM2 "std::array<double, DIM2>":
        double& operator[](size_t)
        double* data()
        size_t size()

    cdef cppclass userdata`'`'DIM "userdata<DIM>":
        double N
        vtype`'`'DIM Pmax
        vtype`'`'DIM sigma
        vtype`'`'DIM alpha
        vtype`'`'DIM beta[DIM]
        double tmax
        double tmin
        void preprocess()

    cdef cppclass pa`'`'DIM "TIN<DIM>":
        pa`'`'DIM`'`'()
        double getEpsilon()
        void setPrecision(double)
        void setUserData(void*)
        void setUB(size_t, double)
        void setShift(vtype`'`'DIM2)
        void optimize()

        vtype`'`'DIM2 xopt
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
        self._thisptr = new pa`'`'DIM`'`'()
        self._ud = new userdata`'`'DIM`'`'()


        for i in range(DIM):
            self._ud.sigma[i] = sigma

        self._thisptr.setUserData(self._ud)

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
    def setChan(self, double[:] alpha, double[:,::1] beta):
        if not (alpha.shape[0] == DIM and alpha.ndim == 1):
            raise ValueError('alpha must be 1D array with DIM elements')

        if not (beta.shape[0] == DIM and beta.shape[1] == DIM and beta.ndim == 2):
            raise ValueError("beta must be DIM()x`'DIM array")

        
        for i in range(DIM):
            self._ud.alpha[i] = alpha[i]

            for j in range(DIM):
                self._ud.beta[i][j] = beta[i][j]

    def setPmax(self, P):
        for i in range(DIM):
            self._ud.Pmax[i] = P
            self._thisptr.setUB(i, P)
    
    def setPrecision(self, epsilon):
        self._thisptr.setPrecision(epsilon)

    def optimize(self):
        self._ud.preprocess();
        self._thisptr.setUB(DIM, self._ud.tmax);

        cdef vtype`'`'DIM2 s
        for i in range(s.size()-1):
            s[i] = .1

        if self._ud.tmin <= 0.1:
            s[s.size()-1] = .1 - self._ud.tmin

        self._thisptr.setShift(s);
        #self._thisptr.setShift(1);
        self._thisptr.optimize()

    def result(self):
        t = self._thisptr
        return np.array([(t.optval, self.xopt(), t.statusStr, t.iter, t.lastUpdate, t.getEpsilon(), t.runtime, t.max_queue_size, t.data_size)], dtype = self.result_dt)

    def getResultDt(self):
        return self.result_dt

    @wraparound(False)
    @boundscheck(False)
    def xopt(self):
        cdef size_t size = self._thisptr.xopt.size() - 1 # don't return auxiliary variable
        cdef cnp.ndarray xopt = np.empty((size,))

        for i in range(size):
            xopt[i] = self._thisptr.xopt[i]

        return xopt
