function buildTrwsMex_time
%mex command to build trwsMex_time

mex src/trwsMex_time.cpp src/ordering.cpp src/MRFEnergy.cpp src/treeProbabilities.cpp src/minimize.cpp -output trwsMex_time -largeArrayDims
