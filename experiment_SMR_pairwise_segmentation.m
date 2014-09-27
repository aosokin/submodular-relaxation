setupSMR;

%% dataset
if ~exist('dataset_alahariPAMI10.mat', 'file')
    createDataset_alahariPAMI10;
end
load('dataset_alahariPAMI10.mat', 'dataset');

%% Run methods on stereo instances
curObjects = 5 : 10;
for iObject = curObjects
    nRows = dataset{iObject}.nRows;
    nCols = dataset{iObject}.nCols;
    nLabels = dataset{iObject}.nLabels;
    dataCost = dataset{iObject}.unary;
    pottsCost = dataset{iObject}.pairwisePotts;
    
    [ vertCost, horCost, pottsCost ] = separateVertHorCosts( pottsCost, [nRows, nCols] );
    
    defaultOptions = struct;
    defaultOptions.funcGetPrimalLabeling = @(initLabels) icmPottsMex(dataCost, pottsCost, initLabels, 5);
    defaultOptions.maxIter = 500;
    
    %% TRW-S
    options_trws = struct;
    options_trws.maxIter = defaultOptions.maxIter;
    options_trws.funcEps = -1;
    [ ~, ~, ~, funcPlot_trws{iObject}, primalPlot_trws{iObject}, timePlot_trws{iObject}] = ...
        trwsMex_time(dataCost, pottsCost, (ones(nLabels, nLabels) - eye(nLabels)), options_trws);
    
    %% SMR
    nVars = nRows * nCols;
    oracle = @(dualVars) computeSmrDualDynamic_pairwisePotts(dataCost, pottsCost, dualVars);
    
    % Hanso
    options = defaultOptions;
    options.lbfgsVectorUpdateNum = 10;
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    [ ~, ~, ~, ~, timePlot_SMR_Hanso{iObject}, funcPlot_SMR_Hanso{iObject}, primalPlot_SMR_Hanso{iObject}] = ...
        maximizeHanso( nVars, oracle, options );
    
    % LMBM
    options = defaultOptions;
    options.maxBundleSize = 10;
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    [ ~, ~, ~, ~, timePlot_SMR_LMBM{iObject}, funcPlot_SMR_LMBM{iObject}, primalPlot_SMR_LMBM{iObject}] = ...
        maximizeLMBM( nVars, oracle, options );
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    
    % Bundle method - fixed size of bundle
    options = defaultOptions;
    options.maxBundleSize = 100;
    options.gammaStepSize = 0.01;
    options.nullSeriosStepRatio = 0.1;
    options.maxInverseStepSize = 100;
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    [ ~, ~, ~, ~, timePlot_SMR_bundle{iObject}, funcPlot_SMR_bundle{iObject}, primalPlot_SMR_bundle{iObject}] = ...
        maximizeBundleMethodFixedSize( nVars, oracle, options );
    
    % Bundle method - Kiwiel's aggregation scheme
    options = defaultOptions;
    options.maxBundleSize = 100;
    options.gammaStepSize = 0.02;
    options.nullSeriosStepRatio = 0.001;
    options.maxInverseStepSize = 1000;
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    [ ~, ~, ~, ~, timePlot_SMR_bundleAggr{iObject}, funcPlot_SMR_bundleAggr{iObject}, primalPlot_SMR_bundleAggr{iObject}] = ...
        maximizeBundleMethodAggregated( nVars, oracle, options );
    
    % Subgradient method
    options = defaultOptions;
    options.gamma = 0.3;
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    [ ~, ~, ~, ~, timePlot_SMR_subgradient{iObject}, funcPlot_SMR_subgradient{iObject}, primalPlot_SMR_subgradient{iObject}] = ...
        maximizeSubgradientAdaptive( nVars, oracle, options );
    
    %% DD TRW
    nVars = nRows * nCols * nLabels;
    oracle = @(dualVars) computeDdtrwDual_pairwisePotts(dataCost, vertCost, horCost, dualVars);
    
        % Hanso
    options = defaultOptions;
    options.lbfgsVectorUpdateNum = 10;
    [ ~, ~, ~, ~, timePlot_DDTRW_Hanso{iObject}, funcPlot_DDTRW_Hanso{iObject}, primalPlot_DDTRW_Hanso{iObject}] = ...
        maximizeHanso( nVars, oracle, options );
    
    % LMBM
    options = defaultOptions;
    options.maxBundleSize = 10;
    [ ~, ~, ~, ~, timePlot_DDTRW_LMBM{iObject}, funcPlot_DDTRW_LMBM{iObject}, primalPlot_DDTRW_LMBM{iObject}] = ...
        maximizeLMBM( nVars, oracle, options );
    
    % Bundle method - fixed size of bundle
    options = defaultOptions;
    options.maxBundleSize = 100;
    options.gammaStepSize = 0.01;
    options.nullSeriosStepRatio = 0.05;
    options.maxInverseStepSize = 2000;
    [ ~, ~, ~, ~, timePlot_DDTRW_bundle{iObject}, funcPlot_DDTRW_bundle{iObject}, primalPlot_DDTRW_bundle{iObject}] = ...
        maximizeBundleMethodFixedSize( nVars, oracle, options );
    
    % Bundle method - Kiwiel's aggregation scheme
    options = defaultOptions;
    options.maxBundleSize = 100;
    options.gammaStepSize = 0.01;
    options.nullSeriosStepRatio = 0.002;
    options.maxInverseStepSize = 500;
    [ ~, ~, ~, ~, timePlot_DDTRW_bundleAggr{iObject}, funcPlot_DDTRW_bundleAggr{iObject}, primalPlot_DDTRW_bundleAggr{iObject}] = ...
        maximizeBundleMethodAggregated( nVars, oracle, options );
    
    % Subgradient method
    options = defaultOptions;
    options.gamma = 0.1;
    [ ~, ~, ~, ~, timePlot_DDTRW_subgradient{iObject}, funcPlot_DDTRW_subgradient{iObject}, primalPlot_DDTRW_subgradient{iObject}] = ...
        maximizeSubgradientAdaptive( nVars, oracle, options );
    
end

computeSmrDualDynamic_pairwisePotts_clearGlobal;

[ticsFuncPlotFinal_SMR_subgradient, funcPlotFinal_SMR_subgradient] = averagePlots( timePlot_SMR_subgradient(curObjects), funcPlot_SMR_subgradient(curObjects), 1 );
[ticsPrimalPlotFinal_SMR_subgradient, primalPlotFinal_SMR_subgradient] = averagePlots( timePlot_SMR_subgradient(curObjects), primalPlot_SMR_subgradient(curObjects), -1 );

[ticsFuncPlotFinal_SMR_bundleAggr, funcPlotFinal_SMR_bundleAggr] = averagePlots( timePlot_SMR_bundleAggr(curObjects), funcPlot_SMR_bundleAggr(curObjects), 1 );
[ticsPrimalPlotFinal_SMR_bundleAggr, primalPlotFinal_SMR_bundleAggr] = averagePlots( timePlot_SMR_bundleAggr(curObjects), primalPlot_SMR_bundleAggr(curObjects), -1 );

[ticsFuncPlotFinal_SMR_bundle, funcPlotFinal_SMR_bundle] = averagePlots( timePlot_SMR_bundle(curObjects), funcPlot_SMR_bundle(curObjects), 1 );
[ticsPrimalPlotFinal_SMR_bundle, primalPlotFinal_SMR_bundle] = averagePlots( timePlot_SMR_bundle(curObjects), primalPlot_SMR_bundle(curObjects), -1 );

[ticsFuncPlotFinal_SMR_LMBM, funcPlotFinal_SMR_LMBM] = averagePlots( timePlot_SMR_LMBM(curObjects), funcPlot_SMR_LMBM(curObjects), 1 );
[ticsPrimalPlotFinal_SMR_LMBM, primalPlotFinal_SMR_LMBM] = averagePlots( timePlot_SMR_LMBM(curObjects), primalPlot_SMR_LMBM(curObjects), -1 );

[ticsFuncPlotFinal_SMR_Hanso, funcPlotFinal_SMR_Hanso] = averagePlots( timePlot_SMR_Hanso(curObjects), funcPlot_SMR_Hanso(curObjects), 1 );
[ticsPrimalPlotFinal_SMR_Hanso, primalPlotFinal_SMR_Hanso] = averagePlots( timePlot_SMR_Hanso(curObjects), primalPlot_SMR_Hanso(curObjects), -1 );

[ticsFuncPlotFinal_DDTRW_subgradient, funcPlotFinal_DDTRW_subgradient] = averagePlots( timePlot_DDTRW_subgradient(curObjects), funcPlot_DDTRW_subgradient(curObjects), 1 );
[ticsPrimalPlotFinal_DDTRW_subgradient, primalPlotFinal_DDTRW_subgradient] = averagePlots( timePlot_DDTRW_subgradient(curObjects), primalPlot_DDTRW_subgradient(curObjects), -1 );

[ticsFuncPlotFinal_DDTRW_bundleAggr, funcPlotFinal_DDTRW_bundleAggr] = averagePlots( timePlot_DDTRW_bundleAggr(curObjects), funcPlot_DDTRW_bundleAggr(curObjects), 1 );
[ticsPrimalPlotFinal_DDTRW_bundleAggr, primalPlotFinal_DDTRW_bundleAggr] = averagePlots( timePlot_DDTRW_bundleAggr(curObjects), primalPlot_DDTRW_bundleAggr(curObjects), -1 );

[ticsFuncPlotFinal_DDTRW_bundle, funcPlotFinal_DDTRW_bundle] = averagePlots( timePlot_DDTRW_bundle(curObjects), funcPlot_DDTRW_bundle(curObjects), 1 );
[ticsPrimalPlotFinal_DDTRW_bundle, primalPlotFinal_DDTRW_bundle] = averagePlots( timePlot_DDTRW_bundle(curObjects), primalPlot_DDTRW_bundle(curObjects), -1 );

[ticsFuncPlotFinal_DDTRW_LMBM, funcPlotFinal_DDTRW_LMBM] = averagePlots( timePlot_DDTRW_LMBM(curObjects), funcPlot_DDTRW_LMBM(curObjects), 1 );
[ticsPrimalPlotFinal_DDTRW_LMBM, primalPlotFinal_DDTRW_LMBM] = averagePlots( timePlot_DDTRW_LMBM(curObjects), primalPlot_DDTRW_LMBM(curObjects), -1 );

[ticsFuncPlotFinal_DDTRW_Hanso, funcPlotFinal_DDTRW_Hanso] = averagePlots( timePlot_DDTRW_Hanso(curObjects), funcPlot_DDTRW_Hanso(curObjects), 1 );
[ticsPrimalPlotFinal_DDTRW_Hanso, primalPlotFinal_DDTRW_Hanso] = averagePlots( timePlot_DDTRW_Hanso(curObjects), primalPlot_DDTRW_Hanso(curObjects), -1 );


[ticsFuncPlotFinal_trws, funcPlotFinal_trws] = averagePlots( timePlot_trws(curObjects), funcPlot_trws(curObjects), 1 );
[ticsPrimalPlotFinal_trws, primalPlotFinal_trws] = averagePlots( timePlot_trws(curObjects), primalPlot_trws(curObjects), -1 );

save('resultsSegmentation.mat');


% plot(ticsFuncPlotFinal_trws, funcPlotFinal_trws, 'k', ticsPrimalPlotFinal_trws, primalPlotFinal_trws, 'k--', ...
%      ticsFuncPlotFinal_SMR_Hanso, funcPlotFinal_SMR_Hanso, 'r', ticsPrimalPlotFinal_SMR_Hanso, primalPlotFinal_SMR_Hanso, 'r--', ...
%      ticsFuncPlotFinal_SMR_LMBM, funcPlotFinal_SMR_LMBM, 'b', ticsPrimalPlotFinal_SMR_LMBM, primalPlotFinal_SMR_LMBM, 'b--', ...
%      ticsFuncPlotFinal_SMR_bundle, funcPlotFinal_SMR_bundle, 'g', ticsPrimalPlotFinal_SMR_bundle, primalPlotFinal_SMR_bundle, 'g--', ...
%      ticsFuncPlotFinal_SMR_bundleAggr, funcPlotFinal_SMR_bundleAggr, 'm', ticsPrimalPlotFinal_SMR_bundleAggr, primalPlotFinal_SMR_bundleAggr, 'm--', ...
%      ticsFuncPlotFinal_SMR_subgradient, funcPlotFinal_SMR_subgradient, 'c', ticsPrimalPlotFinal_SMR_subgradient, primalPlotFinal_SMR_subgradient, 'c--');
% 
