function experiment_SMR_pairwise
% experiment_SMR_pairwise conducts the experiment section 4.1 of the journal paper, results: fig. 4
% CAUTION: rerunning the experiment requires a lot of time and a lot of RAM
% If you just want to examine the plots run makePlots_SMR_segmentation.m and makePlots_SMR_stereo.m
%
% Anton Osokin (firstname.lastname@gmail.com),  14.10.2014

setup_SMR;

%% dataset
if ~exist('dataset_alahariPAMI10.mat', 'file')
    createDataset_alahariPAMI10;
end
load('dataset_alahariPAMI10.mat', 'dataset');
resultFileName = 'results_SMR_recomputed.mat';

%% Run methods on segmentation and stereo instances
curObjects = [1 : 6, 8 : 10];
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
    [ ~, ~, ~, dualPlot_trws{iObject}, primalPlot_trws{iObject}, timePlot_trws{iObject}] = ...
        trwsMex_time(dataCost, pottsCost, (ones(nLabels, nLabels) - eye(nLabels)), options_trws);
    
    %% SMR
    nVars = nRows * nCols;
    oracle = @(dualVars) computeSmrDualDynamic_pairwisePotts(dataCost, pottsCost, dualVars);
    
    % Hanso
    options = defaultOptions;
    options.lbfgsVectorUpdateNum = 10;
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    [ ~, ~, ~, ~, timePlot_SMR_Hanso{iObject}, dualPlot_SMR_Hanso{iObject}, primalPlot_SMR_Hanso{iObject}] = ...
        maximizeHanso( nVars, oracle, options );
    
    % LMBM
    options = defaultOptions;
    options.maxBundleSize = 10;
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    [ ~, ~, ~, ~, timePlot_SMR_LMBM{iObject}, dualPlot_SMR_LMBM{iObject}, primalPlot_SMR_LMBM{iObject}] = ...
        maximizeLMBM( nVars, oracle, options );
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    
    % Bundle method - fixed size of bundle
    options = defaultOptions;
    options.maxBundleSize = 100;
    options.gammaStepSize = 0.01;
    options.nullSeriosStepRatio = 0.1;
    options.maxInverseStepSize = 100;
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    [ ~, ~, ~, ~, timePlot_SMR_bundle{iObject}, dualPlot_SMR_bundle{iObject}, primalPlot_SMR_bundle{iObject}] = ...
        maximizeBundleMethodFixedSize( nVars, oracle, options );
    
    % Bundle method - Kiwiel's aggregation scheme
    options = defaultOptions;
    options.maxBundleSize = 100;
    options.gammaStepSize = 0.02;
    options.nullSeriosStepRatio = 0.001;
    options.maxInverseStepSize = 1000;
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    [ ~, ~, ~, ~, timePlot_SMR_bundleAggr{iObject}, dualPlot_SMR_bundleAggr{iObject}, primalPlot_SMR_bundleAggr{iObject}] = ...
        maximizeBundleMethodAggregated( nVars, oracle, options );
    
    % Subgradient method
    options = defaultOptions;
    options.gamma = 0.3;
    computeSmrDualDynamic_pairwisePotts_clearGlobal;
    [ ~, ~, ~, ~, timePlot_SMR_subgradient{iObject}, dualPlot_SMR_subgradient{iObject}, primalPlot_SMR_subgradient{iObject}] = ...
        maximizeSubgradientAdaptive( nVars, oracle, options );
    
    %% DD TRW
    nVars = nRows * nCols * nLabels;
    oracle = @(dualVars) computeDdtrwDual_pairwisePotts(dataCost, vertCost, horCost, dualVars);
    
        % Hanso
    options = defaultOptions;
    options.lbfgsVectorUpdateNum = 10;
    [ ~, ~, ~, ~, timePlot_DDTRW_Hanso{iObject}, dualPlot_DDTRW_Hanso{iObject}, primalPlot_DDTRW_Hanso{iObject}] = ...
        maximizeHanso( nVars, oracle, options );
    
    % LMBM
    options = defaultOptions;
    options.maxBundleSize = 10;
    [ ~, ~, ~, ~, timePlot_DDTRW_LMBM{iObject}, dualPlot_DDTRW_LMBM{iObject}, primalPlot_DDTRW_LMBM{iObject}] = ...
        maximizeLMBM( nVars, oracle, options );
    
    % Bundle method - fixed size of bundle
    options = defaultOptions;
    options.maxBundleSize = 100;
    options.gammaStepSize = 0.01;
    options.nullSeriosStepRatio = 0.05;
    options.maxInverseStepSize = 2000;
    [ ~, ~, ~, ~, timePlot_DDTRW_bundle{iObject}, dualPlot_DDTRW_bundle{iObject}, primalPlot_DDTRW_bundle{iObject}] = ...
        maximizeBundleMethodFixedSize( nVars, oracle, options );
    
    % Bundle method - Kiwiel's aggregation scheme
    options = defaultOptions;
    options.maxBundleSize = 100;
    options.gammaStepSize = 0.01;
    options.nullSeriosStepRatio = 0.002;
    options.maxInverseStepSize = 500;
    [ ~, ~, ~, ~, timePlot_DDTRW_bundleAggr{iObject}, dualPlot_DDTRW_bundleAggr{iObject}, primalPlot_DDTRW_bundleAggr{iObject}] = ...
        maximizeBundleMethodAggregated( nVars, oracle, options );
    
    % Subgradient method
    options = defaultOptions;
    options.gamma = 0.1;
    [ ~, ~, ~, ~, timePlot_DDTRW_subgradient{iObject}, dualPlot_DDTRW_subgradient{iObject}, primalPlot_DDTRW_subgradient{iObject}] = ...
        maximizeSubgradientAdaptive( nVars, oracle, options );

    
    save(resultFileName, 'timePlot_trws' ,'dualPlot_trws', 'primalPlot_trws', ...
               'timePlot_SMR_Hanso', 'dualPlot_SMR_Hanso', 'primalPlot_SMR_Hanso', ...
               'timePlot_SMR_LMBM', 'dualPlot_SMR_LMBM', 'primalPlot_SMR_LMBM', ...
               'timePlot_SMR_bundle', 'dualPlot_SMR_bundle', 'primalPlot_SMR_bundle', ...
               'timePlot_SMR_bundleAggr', 'dualPlot_SMR_bundleAggr', 'primalPlot_SMR_bundleAggr', ...
               'timePlot_SMR_subgradient', 'dualPlot_SMR_subgradient', 'primalPlot_SMR_subgradient', ...
               'timePlot_DDTRW_Hanso', 'dualPlot_DDTRW_Hanso', 'primalPlot_DDTRW_Hanso', ...
               'timePlot_DDTRW_LMBM', 'dualPlot_DDTRW_LMBM', 'primalPlot_DDTRW_LMBM', ...
               'timePlot_DDTRW_bundle', 'dualPlot_DDTRW_bundle', 'primalPlot_DDTRW_bundle', ...
               'timePlot_DDTRW_bundleAggr', 'dualPlot_DDTRW_bundleAggr', 'primalPlot_DDTRW_bundleAggr', ...
               'timePlot_DDTRW_subgradient', 'dualPlot_DDTRW_subgradient', 'primalPlot_DDTRW_subgradient');
end
computeSmrDualDynamic_pairwisePotts_clearGlobal;

%% plot the results
makePlots_SMR_segmentation(resultFileName);
makePlots_SMR_stereo(resultFileName);

end
