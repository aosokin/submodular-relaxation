This software implements several methods for MRF energy minimization based on the Lagrangian relaxation approach. 
Each method consists of two parts: the convex optimization routine to maximize the dual and the oracle to compute the value and the subgradient of the dual given some value of the dual variables (Lagrange multipliers).

We provide the following oracles:
1) SMR oracle applicable for pairwise associative MRFs;
2) SMR oracle speeded up by the usage of dynamic graph cuts (pairwise associative MRFs);
3) DD TRW oracle applicable for pairwise associative MRFs [27];
4) SMR oracle applicable for energies with high-order robust P^n Potts potentials;
5) DD TRW oracle (CWD) applicable for energies with high-order robust P^n Potts potentials;
6) NSMR oracle applicable for pairwise non-associative MRFs;
7) SMR oracle applicable for pairwise non-associative MRFs via the "subtraction trick".

We provide the following optimization routines:
1) subgradient method with adaptive stepsize [14];
2) bundle method with fixed size of the bundle [14];
3) bundle method with the aggregated bundle (Kiwiel's rule) [14, 16];
4) LMBM method [10];
5) L-BFGS for non-smooth function (Hanso library [30]).

If you use this code, please, cite the following paper in any resulting publication:

Anton Osokin, Dmitry Vetrov, Submodular Relaxation for Inference in Markov Random Fields,
IEEE Transactions on Pattern Analysis and Machine Intelligence, 37(7): 1347-1359, July 2015.

The full text is available on arXiv.org (with some typos corrected):
http://arxiv.org/abs/1501.03771

Anton Osokin, (firstname.lastname@gmail.com)

INSTALLATION
-----------------------------
1) Download the package
2) Run ./setup_SMR from the package folder (requires internet connection to download the GCO library)
3) Run ./example_SMR 

The code was tested under 
- Win7-x64 using MATLAB R2014a and MSVC 2012 + Intel Visual Fortran Composer XE 2013;
- ubuntu-12.04-x64 using MATLAB R2012a and gcc-4.4

Fortran compilers are required exclusively for the LMBM optimization method.
You can use all the other parts of the package without Fortran compilers.

CONTENT OF THE PACKAGE
-----------------------------
./setup_SMR - script to add all the parts of the package to PATH
./build_SMR - script to build all the mexified functions
./example_SMR - an example of usage
./optimizationMethods - code for the included optimization methods
./oracles - code for the oracle functions
./experiments - scripts to reproduce some experiments from the journal paper
./data - functions to work with the datasets used in the journal paper
./mexWrappers - mexified functions including third-party software
./utils - some helper functions

LICENSE
-----------------------------
This software is distributed under MIT license that allows free usage of this code for whatever purposes.

Please note, that some third-party software included in this package are distributed under RESEARCH ONLY license. Please contact the authors of those packages if you interested in using the software for non-research purposes.

THIRD-PARTY SOFTWARE
-----------------------------
This package includes a number of third-party libraries. Please cite the corresponding papers if you are using any of those.

1) Boykov-Kolmogorov max-flow/min-cut algorithm:
Y. Boykov and V. Kolmogorov, An experimental comparison of Min-Cut/Max-Flow algorithms for energy minimization in vision, 
IEEE TPAMI, 26(9):1124-1137, 2004.

2) A dynamic version of Boykov-Kolmogorov max-flow/min-cut algorithm:
Y. Boykov and V. Kolmogorov, An experimental comparison of Min-Cut/Max-Flow algorithms for energy minimization in vision, 
IEEE TPAMI, 26(9):1124-1137, 2004.
P. Kohli and P. Torr, Dynamic graph cuts for efficient inference in markov random fields, IEEE TPAMI, 29(12):2079-2088, 2007.

3) QPBO algorithm for energy minimization:
V. Kolmogorov and C. Rother. Minimizing non-submodular functions with graph cuts - a review. IEEE TPAMI, 29(7):1274-1279, July 2007

4) TRW-S algorithm for energy minimization:
V. Kolmogorov, Convergent tree-reweighted message passing for energy minimization, IEEE TPAMI, vol. 28, no. 10, pp. 1568–1583, 2006.

5) Alpha-expansion algorithm for energy minimization 
Y. Boykov, O. Veksler, R.Zabih, Efficient Approximate Energy Minimization via Graph Cuts, IEEE TPAMI, 20(12):1222-1239, 2001.

6) Hanso - an implementation of L-BFGS method for non-smooth functions. 
A. S. Lewis and M. L. Overton, Nonsmooth optimization via quasi-newton methods, Mathematical Programming, vol. 141, no. 1-2, pp. 135–163, 2013.

7) LMBM optimization method.
N. Haarala, K. Miettinen, and M. M. Makela, Globally convergent limited memory bundle method for large-scale nonsmooth optimization, Mathematical Programming, vol. 109, no. 1, pp. 181–205, 2007.

8) Alpha-expansion for robust P^n Potts potentials.
P. Kohli, L. Ladicky, P. H. S. Torr. 
Robust Higher Order Potentials for Enforcing Label Consistency
International Journal of Computer Vision, 82(3):302-324, 2009.




