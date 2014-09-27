function [ bestPoint, bestVal, bestPrimal, bestLabeling, timePlot, funcPlot, primalPlot ] = maximizeHanso( numVars, oracle, options )
% maximizeHanso maximizes the dual function L-BFGS method form Hanso library:
%   http://www.cs.nyu.edu/overton/software/hanso/hanso2_02.zip
%
% [ bestPoint, bestVal, timePlot, funcPlot, primalPlot ] = maximizeHanso( oracle, options );
%
% INPUT:
%
%
% OUTPUT:
%
%
% Anton Osokin (firstname.lastname@gmail.com),  16.05.2013

if ~exist( 'options', 'var');
    options = struct;
end
options = setDefaultField(options, 'maxIter', 100);
options = setDefaultField(options, 'maxOracleCalls', inf);
options = setDefaultField(options, 'maxTime', inf);
options = setDefaultField(options, 'funcGetPrimalLabeling', @(x) nan);
options = setDefaultField(options, 'initialPoint', []);
options = setDefaultField(options, 'verbose', 'iter');
options = setDefaultField(options, 'lbfgsVectorUpdateNum', []); 



global globalBestPoint
globalBestPoint = nan(numVars, 1);
global globalBestVal;
globalBestVal = -inf;
global globalBestPrimal
globalBestPrimal = inf;
global globalBestLabeling
globalBestLabeling = [];

global globalTimePlot
globalTimePlot = nan(options.maxIter, 1);
global globalFuncPlot
globalFuncPlot = nan(options.maxIter, 1);
global globalPrimalPlot
globalPrimalPlot = nan(options.maxIter, 1);

global globalTimeStart
globalTimeStart = tic;
global globalTimeGarbage
globalTimeGarbage = 0;
global globalNumOracleCalls
globalNumOracleCalls = 0;


if isempty(options.initialPoint)
    initialPoint = -(1e-3)*ones(numVars, 1);
else
    initialPoint = options.initialPoint;
end

hanso_pars = struct;
hanso_pars.nvar = numVars;
hanso_pars.fgname = @(x) inverseOracle(x, oracle, options.funcGetPrimalLabeling);



hanso_options = struct;
hanso_options.x0 = initialPoint;
hanso_options.maxit = options.maxIter;
hanso_options.samprad = [];
hanso_options.cpumax = options.maxTime;
if ~isempty(options.lbfgsVectorUpdateNum)
    hanso_options.nvec = options.lbfgsVectorUpdateNum;
end
if strcmpi(options.verbose, 'iter')
    hanso_options.prtlevel = 2;
elseif strcmpi(options.verbose, 'final')
    hanso_options.prtlevel = 1;
elseif strcmpi(options.verbose, 'none')
    hanso_options.prtlevel = 0;
else
    warning('maximizeHanso:wrongOptionsVerbose', 'Unknown value of options.verbose, expected iter, final, or none');
end

[hansoPoint, hansoVal] = hanso(hanso_pars, hanso_options);

% bestVal, timePlot, funcPlot, primalPlot
bestPoint = hansoPoint;
bestVal = -hansoVal;
if abs( globalBestVal - bestVal ) > 1e-5 * (max(abs(globalBestVal), abs(bestVal)) + 1e-5)
    warning('maximizeHanso:maxLowerBoundMismatch', ['Best LBs computed by Hanso and explicitly do not match, diff: ', num2str( globalBestVal - bestVal )]);
end

timePlot = globalTimePlot;
funcPlot = globalFuncPlot;
primalPlot = globalPrimalPlot;

bestPrimal = globalBestPrimal;
bestLabeling = globalBestLabeling;

clear global globalBestPoint
clear global globalBestVal
clear global globalBestPrimal
clear global globalBestLabeling

clear global globalTimePlot
clear global globalFuncPlot
clear global globalPrimalPlot

clear global globalTimeStart
clear global globalTimeGarbage
clear global globalNumOracleCalls

end

function s = setDefaultField(s, f, value)
    if ~isfield(s, f)
        s.(f) = value;
    end
end

function [fValNegative, subgradientNegative] = inverseOracle(x, oracle, getPrimalLabeling)
    [fVal, subgradient, primalEstimate] = oracle(x);
    fValNegative = -fVal;
    subgradientNegative = -subgradient;
    
    curStart = tic;
    
    global globalBestPoint
    global globalBestVal
    global globalBestPrimal
    global globalBestLabeling

    global globalTimePlot
    global globalFuncPlot
    global globalPrimalPlot

    global globalTimeStart
    global globalTimeGarbage
    global globalNumOracleCalls

    globalNumOracleCalls = globalNumOracleCalls + 1;
    globalFuncPlot( globalNumOracleCalls ) = fVal;
    if fVal > globalBestVal
        globalBestVal = fVal;
        globalBestPoint = x;
    end
    
    [curPrimal, newLabeling] = getPrimalLabeling( primalEstimate );
    if curPrimal < globalBestPrimal
        globalBestPrimal = curPrimal;
        globalBestLabeling = newLabeling;
    end
    globalPrimalPlot( globalNumOracleCalls ) = curPrimal;
    
    globalTimePlot( globalNumOracleCalls ) = nan;
    
    curGarbage = toc( curStart );
    globalTimeGarbage = globalTimeGarbage + curGarbage;
    
    globalTimePlot( globalNumOracleCalls ) = toc( globalTimeStart ) - globalTimeGarbage;
end





