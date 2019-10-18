Mixed Monotonic Programming for Fast Global Optimization
==================

This code package is related to the following scientific article:

Bho Matthiesen, Christoph Hellings, Eduard A. Jorswieck, and Wolfgang Utschick "[Mixed Monotonic Programming for Fast Global Optimization](http://arxiv.org/abs/1910.07853)," submitted to IEEE Transactions on Signal Processing.


## Contents

This is the main directory for the source code. Algorithms and support code is located in `lib/`, application examples that were used numerical example in the paper are in `appl/`.

Building all example applications except python interfaces can be achieved by running `make` in this directory. To build only parts of the source codeplease run make in the respective directory. In this case, please make sure to first compile the code in `lib/`.

## Requirements

This code was compiled and tested with GNU Make 3.82 and GCC 8.2.0. Some application examples require [Gurobi](http://www.gurobi.com/), [Mosek](https://www.mosek.com/), and [Intel MKL](https://software.intel.com/mkl) which are available free of charge for reserach in academic institutions. The simulations for the article cited above were done with Gurobi 8.0.1, Mosek 8.1.0.81, and MKL 2019.

## Acknowledgements

This research was supported in part by the Deutsche Forschungsgemeinschaft (DFG) in the [Collaborative Research Center 912 "Highly Adaptive Energy-Efficient Computing"](https://tu-dresden.de/ing/forschung/sfb912).

We thank the Center for Information Services and High Performance Computing (ZIH) at TU Dresden for generous allocations of computer time.


## License and Referencing

This program is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

