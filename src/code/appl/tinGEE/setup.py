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


from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

exts = [Extension(
            name='geePy',
            sources = ['geePy.pyx'],
            include_dirs = ['.','../../lib',numpy.get_include()],
            extra_objects = ['../../lib/util.o'],
            libraries = ['m',],
            language = "c++",
            extra_compile_args = ["-std=c++17", '-fpic', "-O3", '-march=haswell', '-mtune=haswell', '-DNDEBUG', '-m64'],
            extra_link_args = ["-std=c++17"]),
         Extension(
            name='dinkelbachPy',
            sources = ['dinkelbachPy.pyx'],
            include_dirs = ['.','../../lib',numpy.get_include()],
            extra_objects = ['../../lib/util.o'],
            libraries = ['m',],
            language = "c++",
            extra_compile_args = ["-std=c++17", '-fpic', "-O3", '-march=haswell', '-mtune=haswell', '-DNDEBUG', '-m64'],
            extra_link_args = ["-std=c++17"]),
        ]

setup(
    ext_modules = cythonize(exts)
)
