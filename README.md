Mixed Monotonic Programming for Fast Global Optimization
==================

This code package is related to the following scientific article:

Bho Matthiesen, Christoph Hellings, Eduard A. Jorswieck, and Wolfgang Utschick "[Mixed Monotonic Programming for Fast Global Optimization](http://arxiv.org/abs/1910.07853)," submitted to IEEE Transactions on Signal Processing.


## Abstract of the Article

While globally optimal solutions to convex programs can be computed efficiently in polynomial time, this is, in general, not possible for nonconvex optimization problems. Therefore, locally optimal approaches or other efficient suboptimal heuristics are usually applied for practical implementations. However, there is also a strong interest in computing globally optimal solutions of nonconvex problems in offline simulations in order to benchmark the faster suboptimal algorithms. Global solutions often rely on monotonicity properties. A common approach is to reformulate problems into a canonical form of a monotonic optimization problem, where the monotonicity becomes evident, but this often comes at the cost of nested optimizations, increased numbers of variables, and/or slow convergence. The framework of mixed monotonic programming (MMP) proposed in this paper is a more direct approach that exploits hidden monotonicity properties without performance-deteriorating reformulations. By means of a wide range of application examples from the area of signal processing for communications (including energy efficiency for green communications, resource allocation in interference networks, scheduling for fairness and quality of service, as well as beamformer design in multiantenna systems), we demonstrate that the novel MMP approach leads to tremendous complexity reductions compared to state-of-the-art methods for global optimization.

## Requirements & Contents of the Code Package

This code was compiled and tested with GNU Make 3.82 and GCC 8.2.0. Some application examples require [Gurobi](http://www.gurobi.com/), [Mosek](https://www.mosek.com/), and [Intel MKL](https://software.intel.com/mkl) which are available free of charge for reserach in academic institutions. The simulations for the article cited above were done with Gurobi 8.0.1, Mosek 8.1.0.81, and MKL 2019.

The complete source code is available in `src/code/`. Prior to compilation, please update the variables in the Makefile according to your needs. The simulations were conducted on TU Dresden's HPC. Sample [SLURM](https://www.schedmd.com/) `sbatch` files that illustrate the usage of the code are available in `src/slurm/` and some other run scripts in `src/noslurm`. All input data is stored in `data` and raw results plus evaluation scripts are in `results/`.

## Acknowledgements

This research was supported in part by the Deutsche Forschungsgemeinschaft (DFG) in the [Collaborative Research Center 912 "Highly Adaptive Energy-Efficient Computing"](https://tu-dresden.de/ing/forschung/sfb912).

We thank the Center for Information Services and High Performance Computing (ZIH) at TU Dresden for generous allocations of computer time.


## License and Referencing

This program is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

