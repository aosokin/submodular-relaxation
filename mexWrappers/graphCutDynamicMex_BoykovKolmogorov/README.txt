This software implements the MATLAB wrapper for Boykov-Kolmogorov max-flow/min-cut algorithm.

Anton Osokin, (firstname.lastname@gmail.com)
24.09.2014

Please cite the following paper in any resulting publication:

Y. Boykov and V. Kolmogorov, An experimental comparison of Min-Cut/Max-Flow algorithms for energy minimization in vision, 
IEEE TPAMI, 26(9):1124–1137, 2004.

P. Kohli and P. Torr, Dynamic graph cuts for efficient
inference in markov random fields,
IEEE TPAMI, 29(12):2079–2088, 2007.

PACKAGE
-----------------------------

./graphCutDynamicMex.cpp, ./updateGraphCutDynamicMex.cpp, ./deletegraphCutDynamicMex.cpp, ./graphCutMemory.h, ./graphCutMemory.cpp, , ./graphCutMex.h,  - the C++ code of the wrapper

./build_graphCutDynamicMex.m - function to build the wrapper

./graphCutDynamicMex.m, ./updateGraphCutDynamicMex.m, ./deleteGraphCutDynamicMex.m - the description of the implemented functions

./example_graphCutDynamicMex.m - the example of usage

./maxflow-v3.03.src - C++ code by Vladimir Kolmogorov (the code was slightly modified)
http://pub.ist.ac.at/~vnk/software/maxflow-v3.03.src.zip

./graphCutDynamicMex.mexw64, ./updateGraphCutDynamicMex.mexw64, ./deleteGraphCutDynamicMex.mexw64 - binary files for the MEX-functions compiled using MATLAB R2014a + MSVC 2012

USING THE CODE
-----------------------------

0) Install MATLAB and one of the supported compilers

1) Run build_graphCutDynamicMex.m

2) Run example_graphCutDynamicMex.m to check if the code works

The code was tested using MSVC 2012 and MATLAB 2014a

