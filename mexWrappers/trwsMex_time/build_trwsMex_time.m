function build_trwsMex_time
% build_trwsMex_time builds package trwsMex_time
%
% Anton Osokin (firstname.lastname@gmail.com), 24.09.2014

mex src/trwsMex_time.cpp src/ordering.cpp src/MRFEnergy.cpp src/treeProbabilities.cpp src/minimize.cpp -output trwsMex_time -largeArrayDims
