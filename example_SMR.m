% example of usage of the SMR package

setup_SMR;
load('example_SMR.mat', 'exampleEnergy');

nRows = exampleEnergy.nRows;
nCols = exampleEnergy.nCols;
nLabels = exampleEnergy.nLabels;
dataCost = double(exampleEnergy.unary);
pottsCost = exampleEnergy.pairwisePotts;

options = struct;
options.funcGetPrimalLabeling = @(initLabels) icmPottsMex(dataCost, pottsCost, initLabels, 5);
options.verbose = 'iter';
options.maxIter = 100;

nVars = nRows * nCols;
oracle = @(dualVars) computeSmrDualDynamic_pairwisePotts(dataCost, pottsCost, dualVars);

[ bestPoint, bestVal, bestPrimal, bestLabeling, timePlot, dualPlot, primalPlot] = ...
    maximizeHanso( nVars, oracle, options );
computeSmrDualDynamic_pairwisePotts_clearGlobal;

figure(1);
plot(timePlot(10:end), dualPlot(10:end), 'b', timePlot(10:end), primalPlot(10:end), 'r');
ylabel('energy and lower bound');
xlabel('time (seconds)');
title('SMR energies and lower bounds');
legend({'lower bound', 'energy'});

figure(2);
imshow( reshape(bestLabeling, [nRows, nCols]), [] );
