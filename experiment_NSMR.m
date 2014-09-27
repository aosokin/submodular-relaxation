setupSMR;

%% dataset
createDataset_NSMR;
load('dataset_NSMR.mat', 'dataset');

%%  run TRW-S, NSMR, SMR + subtraction, DD TRW
funcPlot_trws = cell(length(dataset), 1);
funcPlot_ddtrw = cell(length(dataset), 1);
funcPlot_nsmr = cell(length(dataset), 1);
funcPlot_subtraction = cell(length(dataset), 1);

timePlot_trws = cell(length(dataset), 1);
timePlot_ddtrw = cell(length(dataset), 1);
timePlot_nsmr = cell(length(dataset), 1);
timePlot_subtraction = cell(length(dataset), 1);

primalPlot_trws = cell(length(dataset), 1);
primalPlot_ddtrw = cell(length(dataset), 1);
primalPlot_nsmr = cell(length(dataset), 1);
primalPlot_subtraction = cell(length(dataset), 1);

for iObject = 1 : length(dataset)
    nRows = dataset{iObject}.nRows;
    nCols = dataset{iObject}.nCols;
    nVars = nRows * nCols;
    nLabels = dataset{iObject}.nLabels;
    dataCost = dataset{iObject}.unary;
    pottsCost = dataset{iObject}.pairwisePotts;
    
    
    options = struct;
    options.funcGetPrimalLabeling = @(initLabels) icmPottsMex(dataCost, pottsCost, initLabels, 5);
    options.maxIter = 200;
    
    % SMR with subtraction
    oracle_subtraction = @(dualVars) computeSmrDualSubtraction_pairwisePotts(dataCost, pottsCost, dualVars);
    [ bestPoint_subtraction, bestVal_subtraction, bestPrimal_subtraction, bestLabeling_subtraction, timePlot_subtraction{iObject}, funcPlot_subtraction{iObject}, primalPlot_subtraction{iObject} ] = ...
        maximizeHanso( nVars, oracle_subtraction, options );
    
    % NSMR
    oracle_nsmr = @(dualVars) computeNsmrDual_pairwisePotts(dataCost, pottsCost, dualVars);
    [ bestPoint_nsmr, bestVal_nsmr, bestPrimal_nsmr, bestLabeling_nsmr, timePlot_nsmr{iObject}, funcPlot_nsmr{iObject}, primalPlot_nsmr{iObject} ] = ...
        maximizeHanso( nVars, oracle_nsmr, options );
    
    % TRW-S
    options_trws = struct;
    options_trws.maxIter = options.maxIter;
    options_trws.funcEps = -1;
    [ bestLabeling_trws,  bestPrimal_trws, bestVal_trws, funcPlot_trws{iObject}, primalPlot_trws{iObject}, timePlot_trws{iObject}] = ...
        trwsMex_time(dataCost, pottsCost, (ones(nLabels, nLabels) - eye(nLabels)), options_trws);
    
    % DD TRW
    [ vertCost, horCost ] = separateVertHorCosts( pottsCost, [nRows; nCols] );
    oracle_ddtrw = @(dualVars) computeDdtrwDualPotts(dataCost, vertCost, horCost, dualVars);
    [ bestPoint_ddtrw, bestVal_ddtrw, bestPrimal_ddtrw, bestLabeling_ddtrw, timePlot_ddtrw{iObject}, funcPlot_ddtrw{iObject}, primalPlot_ddtrw{iObject} ] = ...
        maximizeHanso( nVars * nLabels, oracle_ddtrw, options );
end

for iObject = 1 : length(dataset)
    minPrimal = min( primalPlot_trws{iObject}(:) );
    minPrimal = min( minPrimal, min( primalPlot_nsmr{iObject}(:)) );
    minPrimal = min( minPrimal, primalPlot_subtraction{iObject} );
    minPrimal = min( minPrimal, min( primalPlot_ddtrw{iObject}(:)) );
    
    gapTrws(iObject) = minPrimal - max(funcPlot_trws{iObject}(:));
    gapNsmr(iObject) = minPrimal - max(funcPlot_nsmr{iObject}(:));
    gapSubtraction(iObject) = minPrimal - max(funcPlot_subtraction{iObject});
    gapDdtrw(iObject) = minPrimal - max(funcPlot_ddtrw{iObject}(:));
end

[medianTrws, lowerQuartileTrws, upperQuatileTrws] = getQuartile(gapTrws);
[medianNsmr, lowerQuartileNsmr, upperQuatileNsmr] = getQuartile(gapNsmr);
[medianSubtraction, lowerQuartileSubtraction, upperQuatileSubtraction] = getQuartile(gapSubtraction);
[medianDdtrw, lowerQuartileDdtrw, upperQuatileDdtrw] = getQuartile(gapDdtrw);







