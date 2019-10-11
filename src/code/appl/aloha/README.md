Mixed Monotonic Programming for Fast Global Optimization
==================

This code package is related to the following scientific article:

Bho Matthiesen, Christoph Hellings, Eduard A. Jorswieck, and Wolfgang Utschick "Mixed Monotonic Programming for Fast Global Optimization," submitted to IEEE Transactions on Signal Processing.


## Contents

Probability optimization for slotted ALOHA. Corresponds to Section IV-F.

* `aloha.{cpp,h}`, `aloha4.cpp`: MMP Aloha.
* `aloha_double.cpp`: Polyblock Aloha from [1].
* `aloha_double_dm4.cpp`, `aloha_double_dm.{cpp,h}`: Branch-and-bound implementation of Aloha for DM bounds.
* `aloha_gp.{cpp,h}`, `aloha_gp4.cpp`: Geometric programming based solution from [1].
* `aloha_helper.h`: Implementation helper for applications and HDF interface to read input data and store results.
* `dgopt.{c,h}`: Dual GP solver from Mosek distribution (minor modifications).
* `Mosek.{cpp,h}`: Mosek helper functions.

* `collect_aloha.py`: Collect results in a single file.
* `aloha_createData_{3,4}.py`: Generate input data for numerical expermient.

References:

1. Y. J. Zhang, L. P. Qian, and J. Huang, Monotonic Optimization in Communication and Networking Systems, ser. FnT Netw. Boston, MA, USA: Now, 2012, vol. 7, no. 1, ch. 7.

## Requirements

This code was compiled and tested with GNU Make 3.82 and GCC 8.2.0. This application example requires [Mosek](https://www.mosek.com/), [Intel MKL](https://software.intel.com/mkl), and [HDF5](https://www.hdfgroup.org/) which are available free of charge for reserach in academic institutions. The simulations for the article cited above were done with Mosek 8.1.0.81 and MKL 2019.

## Acknowledgements

This research was supported in part by the Deutsche Forschungsgemeinschaft (DFG) in the [Collaborative Research Center 912 "Highly Adaptive Energy-Efficient Computing"](https://tu-dresden.de/ing/forschung/sfb912).

We thank the Center for Information Services and High Performance Computing (ZIH) at TU Dresden for generous allocations of computer time.


## License and Referencing

This program is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

