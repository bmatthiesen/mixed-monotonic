Mixed Monotonic Programming for Fast Global Optimization
==================

This code package is related to the following scientific article:

Bho Matthiesen, Christoph Hellings, Eduard A. Jorswieck, and Wolfgang Utschick "Mixed Monotonic Programming for Fast Global Optimization," submitted to IEEE Transactions on Signal Processing.


## Contents

Proportional fair rate optimization with scheduling. Corresponds to Section IV-C.

* `tin_pf_dual.cpp`: solves the proportional fair scheduling problem via Lagrange duality and MMP
* `tin_asmapel.cpp`: solves the proportional fair scheduling problem approximately via ASMAPEL
* `ASMAPEL_PA.h`: A-S-MAPEL implementation from L. P. Qian and Y. J. Zhang, “S-MAPEL: Monotonic optimization for non-convex joint power control and scheduling problems,” IEEE Trans. Wireless Commun., vol. 9, no. 5, pp. 1708–1719, May 2010.
* `bits`: files included by the header files in this directory

## Requirements

This code was compiled and tested with GNU Make 3.82 and GCC 8.2.0. This application examples require [Gurobi](http://www.gurobi.com/) and [Intel MKL](https://software.intel.com/mkl) which are available free of charge for reserach in academic institutions. The simulations for the article cited above were done with Gurobi 8.0.1 and MKL 2019.

## Acknowledgements

This research was supported in part by the Deutsche Forschungsgemeinschaft (DFG) in the [Collaborative Research Center 912 "Highly Adaptive Energy-Efficient Computing"](https://tu-dresden.de/ing/forschung/sfb912).

We thank the Center for Information Services and High Performance Computing (ZIH) at TU Dresden for generous allocations of computer time.


## License and Referencing

This program is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

