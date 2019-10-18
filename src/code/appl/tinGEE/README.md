Mixed Monotonic Programming for Fast Global Optimization
==================

This code package is related to the following scientific article:

Bho Matthiesen, Christoph Hellings, Eduard A. Jorswieck, and Wolfgang Utschick "[Mixed Monotonic Programming for Fast Global Optimization](http://arxiv.org/abs/1910.07853)," submitted to IEEE Transactions on Signal Processing.


## Contents

Energy efficiency optimization. Corresponds to Section IV-B.

* `gee_tin.h`: MMP implementation.
* `tinEE_dinkelbach.h`: Dinkelbach + branch-and-bound implementation.

* `gee_tin.cpp`: Test `gee_tin.h`.
* `tinEE_dm.cpp`: Test `tinEE_dm.h`.

* `collect.py`: Collect results.
* `run_gee.py`: Run numerical experiments.
* `setup.py`: Called from Makefile to build python code.
* `test_gee.py`: Simple script to test compiled modules.


* `*.pyx.m4`: Macro files to generate cython interface.

## Building

Run `make` to build standalone applications and `make python` for the python interface.

## Requirements

This code was compiled and tested with GNU Make 3.82, GCC 8.2.0 and the m4 macro processor. The python interface requires cython which is, for example, included in the [Anaconda](https://www.anaconda.com/distribution/) distribution.

## Acknowledgements

This research was supported in part by the Deutsche Forschungsgemeinschaft (DFG) in the [Collaborative Research Center 912 "Highly Adaptive Energy-Efficient Computing"](https://tu-dresden.de/ing/forschung/sfb912).

We thank the Center for Information Services and High Performance Computing (ZIH) at TU Dresden for generous allocations of computer time.


## License and Referencing

This program is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

