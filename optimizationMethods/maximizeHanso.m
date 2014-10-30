function [ bestDualPoint, bestDualValue, bestPrimalValue, bestPrimalLabeling, timePlot, dualPlot, primalPlot ] = maximizeHanso( numVars, oracle, options )
% maximizeHanso maximizes the dual function with L-BFGS method form Hanso library:
%   http://www.cs.nyu.edu/overton/software/hanso/hanso2_02.zip
%
% Please cite the following paper if you use this code:
%   A. S. Lewis and M. L. Overton, Nonsmooth optimization via quasi-newton methods,
%   Mathematical Programming, vol. 141, no. 1-2, pp. 135-163, 2013.
%
% Usage:
% [ bestDualPoint, bestDualValue, bestPrimalValue, bestPrimalLabeling, timePlot, dualPlot, primalPlot ] = maximizeHanso( numVars, oracle, options );
%
% INPUT:
%   numVars - number of variables to optimize over (a single number)
%   oracle - function handle to the function to optimize
%       should be of format [fValue, subgradient, primalEstimate] = oracle(X); where
%           X is a vector of dimensionality numVars, 
%           fValue - value of the function, 
%           subgradient - a subgradient of the optimized function, 
%           primalEstimate - the current estimate of primal labeling
%   options (optional) - a structure with the parameter of the method, fields:
%       maxIter - maximum number of iterations (default: 100)
%       maxTime - maximum running time (default: inf)
%       funcGetPrimalLabeling - function to compute the current estimate of the primal solution;
%           should be of format [ energy, labeling ] = funcGetPrimalLabeling( primalEstimate ); where
%               energy - is the energy value;
%               labeling - the current labeling;
%               primalEstimate - the initialization for the labeling coming from the dual
%           (default: nan for all labelings), consider using icmPottsMex
%       initialPoint - the initialization of the dual variables, (default: -(1e-3)*ones(numVars, 1)).
%       verbose - the verbosity level, could be 'iter', 'final', or 'none' (default: 'iter')
%       lbfgsVectorUpdateNum - the number of vectors to save and use in the limited memory updates 
%           (default: 0 if numVars <= 100, otherwise 10);
%       argTol - argument convergence threshhold (default: 1e-4)
%
% OUTPUT:
%       bestDualPoint - the value of the variables that delivers the maximum-found function value
%       bestDualValue - the function value at bestDualPoint
%       bestPrimalValue - the lowest-found energy
%       bestPrimalLabeling - the labeling that delivers bestPrimalValue
%       timePlot, dualPlot, primalPlot - running time, the dual and the primal values for each oracle call 
%
%   Depends on Hanso v-2.02: optimizationMethods/hanso2_02
%
% Anton Osokin (firstname.lastname@gmail.com),  16.05.2013

if ~exist( 'options', 'var');
    options = struct;
end
options = setDefaultField(options, 'maxIter', 100);
options = setDefaultField(options, 'maxTime', inf);
options = setDefaultField(options, 'funcGetPrimalLabeling', @(x) funcGetPrimalLabeling_default(x));
options = setDefaultField(options, 'initialPoint', []);
options = setDefaultField(options, 'verbose', 'iter');
options = setDefaultField(options, 'lbfgsVectorUpdateNum', []); 
options = setDefaultField(options, 'argTol', 1e-4); 

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
hanso_pars.fgname = @(x) negativeOracle(x, oracle, options.funcGetPrimalLabeling);



hanso_options = struct;
hanso_options.x0 = initialPoint;
hanso_options.maxit = options.maxIter;
hanso_options.samprad = [];
hanso_options.cpumax = options.maxTime;
hanso_options.normtol = options.argTol;

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
bestDualPoint = hansoPoint;
bestDualValue = -hansoVal;
if abs( globalBestVal - bestDualValue ) > 1e-5 * (max(abs(globalBestVal), abs(bestDualValue)) + 1e-5)
    warning('maximizeHanso:maxLowerBoundMismatch', ['Best LBs computed by Hanso and explicitly do not match, diff: ', num2str( globalBestVal - bestDualValue )]);
end

lastNonNan = find(~isnan(globalTimePlot), 1, 'last');
timePlot = globalTimePlot(1 : lastNonNan);
dualPlot = globalFuncPlot(1 : lastNonNan);
primalPlot = globalPrimalPlot(1 : lastNonNan);

bestPrimalValue = globalBestPrimal;
bestPrimalLabeling = globalBestLabeling;

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

function [fValNegative, subgradientNegative] = negativeOracle(x, oracle, getPrimalLabeling)
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





