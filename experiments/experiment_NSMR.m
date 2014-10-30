function [medianTrws, lowerQuartileTrws, upperQuartileTrws, ...
    medianNsmr, lowerQuartileNsmr, upperQuartileNsmr, ...
    medianSubtraction, lowerQuartileSubtraction, upperQuartileSubtraction, ...
    medianDdtrw, lowerQuartileDdtrw, upperQuartileDdtrw] = experiment_NSMR
% experiment_NSMR conducts the experiment from section B, supplemetary material, results: table 3
%
% Anton Osokin (firstname.lastname@gmail.com),  14.10.2014

setup_SMR;

%% dataset
createDataset_NSMR;
load('dataset_NSMR.mat', 'dataset');

%%  run TRW-S, NSMR, SMR + subtraction, DD TRW
dualPlot_trws = cell(length(dataset), 1);
dualPlot_ddtrw = cell(length(dataset), 1);
dualPlot_nsmr = cell(length(dataset), 1);
dualPlot_subtraction = cell(length(dataset), 1);

timePlot_trws = cell(length(dataset), 1);
timePlot_ddtrw = cell(length(dataset), 1);
timePlot_nsmr = cell(length(dataset), 1);
timePlot_subtraction = cell(length(dataset), 1);

primalPlot_trws = cell(length(dataset), 1);
primalPlot_ddtrw = cell(length(dataset), 1);
primalPlot_nsmr = cell(length(dataset), 1);
primalPlot_subtraction = cell(length(dataset), 1);

numObjects = length(dataset);
parfor iObject = 1 : numObjects
    fprintf('Instance %d of %d\n', iObject, numObjects);
    nRows = dataset{iObject}.nRows;
    nCols = dataset{iObject}.nCols;
    nVars = nRows * nCols;
    nLabels = dataset{iObject}.nLabels;
    dataCost = dataset{iObject}.unary;
    pottsCost = dataset{iObject}.pairwisePotts;
    
    
    options = struct;
    options.funcGetPrimalLabeling = @(initLabels) icmPottsMex(dataCost, pottsCost, initLabels, 5);
    options.maxIter = 200;
    options.verbose = 'none';
    
    % SMR with subtraction
    oracle_subtraction = @(dualVars) computeSmrDualSubtraction_pairwisePotts(dataCost, pottsCost, dualVars);
    [ ~, ~, ~, ~, timePlot_subtraction{iObject}, dualPlot_subtraction{iObject}, primalPlot_subtraction{iObject} ] = ...
        maximizeHanso( nVars, oracle_subtraction, options );
    
    % NSMR
    oracle_nsmr = @(dualVars) computeNsmrDual_pairwisePotts(dataCost, pottsCost, dualVars);
    [ ~, ~, ~, ~, timePlot_nsmr{iObject}, dualPlot_nsmr{iObject}, primalPlot_nsmr{iObject} ] = ...
        maximizeHanso( nVars, oracle_nsmr, options );
    
    % TRW-S
    options_trws = struct;
    options_trws.maxIter = options.maxIter;
    options_trws.funcEps = -1;
    options_trws.verbosity = 0;
    
    [ ~, ~, ~, dualPlot_trws{iObject}, primalPlot_trws{iObject}, timePlot_trws{iObject}] = ...
        trwsMex_time(dataCost, pottsCost, (ones(nLabels, nLabels) - eye(nLabels)), options_trws);
    
    % DD TRW
    [ vertCost, horCost ] = separateVertHorCosts( pottsCost, [nRows; nCols] );
    oracle_ddtrw = @(dualVars) computeDdtrwDual_pairwisePotts(dataCost, vertCost, horCost, dualVars);
    [ ~, ~, ~, ~, timePlot_ddtrw{iObject}, dualPlot_ddtrw{iObject}, primalPlot_ddtrw{iObject} ] = ...
        maximizeHanso( nVars * nLabels, oracle_ddtrw, options );
end

% estimate the integrality gaps
gapTrws = nan(numObjects, 1);
gapNsmr = nan(numObjects, 1);
gapSubtraction = nan(numObjects, 1);
gapDdtrw = nan(numObjects, 1);
for iObject = 1 : numObjects
    gapTrws(iObject) = min(primalPlot_trws{iObject}(:)) - max(dualPlot_trws{iObject}(:));
    gapNsmr(iObject) = min( primalPlot_nsmr{iObject}(:)) - max(dualPlot_nsmr{iObject}(:));
    gapSubtraction(iObject) = min( primalPlot_subtraction{iObject}(:)) - max(dualPlot_subtraction{iObject});
    gapDdtrw(iObject) = min( primalPlot_ddtrw{iObject}(:)) - max(dualPlot_ddtrw{iObject}(:));

end

% compute the quartiles
[medianTrws, lowerQuartileTrws, upperQuartileTrws] = getQuartile(gapTrws);
[medianNsmr, lowerQuartileNsmr, upperQuartileNsmr] = getQuartile(gapNsmr);
[medianSubtraction, lowerQuartileSubtraction, upperQuartileSubtraction] = getQuartile(gapSubtraction);
[medianDdtrw, lowerQuartileDdtrw, upperQuartileDdtrw] = getQuartile(gapDdtrw);

% print the results
fprintf('NSMR    gap: lower quartile %f, median %f, upper quartile %f\n', lowerQuartileNsmr, medianNsmr, upperQuartileNsmr);
fprintf('SMR+sub gap: lower quartile %f, median %f, upper quartile %f\n', lowerQuartileSubtraction, medianSubtraction, upperQuartileSubtraction);
fprintf('TRW-S   gap: lower quartile %f, median %f, upper quartile %f\n', lowerQuartileTrws, medianTrws, upperQuartileTrws);
fprintf('DD-TRW  gap: lower quartile %f, median %f, upper quartile %f\n', lowerQuartileDdtrw, medianDdtrw, upperQuartileDdtrw);





