Mixed Monotonic Programming for Fast Global Optimization
==================

This code package is related to the following scientific article:

Bho Matthiesen, Christoph Hellings, Eduard A. Jorswieck, and Wolfgang Utschick "[Mixed Monotonic Programming for Fast Global Optimization](http://arxiv.org/abs/1910.07853)," submitted to IEEE Transactions on Signal Processing.


## Contents

Weighted sum rate maximization in a K-user interference channel. Corresponds to Section IV-A.

* `tin.h`: MMP implementation.
* `tin_fifo.h`: MMP implementation with oldest-first selection.
* `tin_dm.cpp`: BRB implementation with DM bound.
* `tin_PA.h`: Polyblock implementation.
* `tin_mapel.h`: MAPEL algorithm [1].
* `tin_ratespace.h`: Ratespace polyblock algorithm [2], [3], [4].

* `tin.cpp`: Test `tin.h`.
* `tin_dm.h`: Test `tin_dm.h`.
* `tin_ratespace.cpp`: Test `tin_ratespace.h`.
* `tin_PA.cpp`: Test `tin_PA.h`.
* `tin_mapel.cpp`: Test `tin_mapel.h`.

* `run_tin.py`, `run_tin2.py`: Run scripts.
* `collect.py`: Collect results.
* `filelock.py`: Simple locking facilities.
* `setup.py`: Called from Makefile to build python code.
* `test_tin.py`: Simple script to test compiled modules.

* `*.pyx.m4`: Macro files to generate cython interface.

References:

1. L. P. Qian, Y. J. Zhang, and J. Huang, “MAPEL: Achieving global optimality for a non-convex wireless power control problem,” IEEE Trans. Wireless Commun., vol. 8, no. 3, pp. 1553–1563, Mar. 2009.
2. L. Liu, R. Zhang, and K.-C. Chua, “Achieving global optimality for weighted sum-rate maximization in the K-user Gaussian interference channel with multiple antennas,” IEEE Trans. Wireless Commun., vol. 11, no. 5, pp. 1933–1945, May 2012.
3. W. Utschick and J. Brehmer, “Monotonic optimization framework for coordinated beamforming in multicell networks,” IEEE Trans. Signal Process., vol. 60, no. 4, pp. 1899–1909, Apr. 2012.
4. J. Brehmer, Utility Maximization in Nonconvex Wireless Systems, ser. Found. Signal Process., Commun., Netw., W. Utschick, H. Boche, and R. Mathar, Eds. Springer-Verlag, 2012, vol. 5.  


## Building

Run `make` to build standalone applications and `make python` for the python interface.

## Requirements

This code was compiled and tested with GNU Make 3.82, GCC 8.2.0 and the m4 macro processor. This application example requires [Intel MKL](https://software.intel.com/mkl) which is available free of charge for reserach in academic institutions. The simulations for the article cited above were done with MKL 2019. The python interface requires cython which is, for example, included in the [Anaconda](https://www.anaconda.com/distribution/) distribution.

## Acknowledgements

This research was supported in part by the Deutsche Forschungsgemeinschaft (DFG) in the [Collaborative Research Center 912 "Highly Adaptive Energy-Efficient Computing"](https://tu-dresden.de/ing/forschung/sfb912).

We thank the Center for Information Services and High Performance Computing (ZIH) at TU Dresden for generous allocations of computer time.


## License and Referencing

This program is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

