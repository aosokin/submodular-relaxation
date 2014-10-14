function [ bestDualPoint, bestDualValue, bestPrimalValue, bestPrimalLabeling, timePlot, dualPlot, primalPlot ] = maximizeSubgradientAdaptive( numVars, oracle, options )
%maximizeSubgradientAdaptive maximizes the dual function using the bundle method 
% with bundle of limited size.
%
% Please cite the following paper if you use this code:
%   J. H. Kappes, B. Savchynskyy, and C. Schnoerr, 
%   A bundle approach to efficient MAP-inference by Lagrangian relaxation,
%   CVPR, 2012.
%
% Usage:
% [ bestDualPoint, bestDualValue, bestPrimalValue, bestPrimalLabeling, timePlot, dualPlot, primalPlot ] = maximizeSubgradientAdaptive( numVars, oracle, options );
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
%       verbose - the verbosity level, could be 'iter' or 'none' (default: 'iter')
%       gamma - the method parameter \gamma; see the journal paper for the details; (default: 0.1)
%       argTolEps - if the current step size is less than argTolEps the optimization terminates; (default: 1e-20)
%       funcEps - if the estimate of the integrality gap is less than funcEps the optimization terminates; (default: 1e-5)
%
% OUTPUT:
%       bestDualPoint - the value of the variables that delivers the maximum-found function value
%       bestDualValue - the function value at bestDualPoint
%       bestPrimalValue - the lowest-found energy
%       bestPrimalLabeling - the labeling that delivers bestPrimalValue
%       timePlot, dualPlot, primalPlot - running time, the dual and the primal values for each oracle call 
%
%   Depends ./optimizationMethods/bundleMethod
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
options = setDefaultField(options, 'gamma', 0.1);
options = setDefaultField(options, 'argTolEps', 1e-20);
options = setDefaultField(options, 'funcEps', 1e-5);



bestDualPoint = nan(numVars, 1);
bestDualValue = -inf;
bestPrimalValue = inf;
bestPrimalLabeling = [];
timePlot = nan(options.maxIter, 1);
dualPlot = nan(options.maxIter, 1);
primalPlot = nan(options.maxIter, 1);

tStart = tic;
timeGarbage = 0;
numOracleCalls = 0;
reachedMaxIter = true;

if isempty(options.initialPoint)
    initialPoint = -(1e-3)*ones(numVars, 1);
else
    initialPoint = options.initialPoint;
end

currentPoint = initialPoint;

for iIteration = 1 : options.maxIter
    
    [fValue, subgradient, primalEstimate] = oracle( currentPoint );
    dualPlot(iIteration) = fValue;
    timePlot(iIteration) = toc(tStart) - timeGarbage;
    numOracleCalls = numOracleCalls + 1;
    % exclude primal estimation from time 
    tTemp = tic;
    [primalPlot( iIteration ), newLabeling ] = options.funcGetPrimalLabeling(primalEstimate);
%     timeGarbage = timeGarbage + toc(tTemp);
    
    if fValue > bestDualValue
        bestDualValue = fValue;
        bestDualPoint = currentPoint;
    end
    if primalPlot( iIteration ) < bestPrimalValue
        bestPrimalValue = primalPlot( iIteration );
        bestPrimalLabeling = newLabeling;
    end
    
    
   alpha = options.gamma * (bestPrimalValue - fValue) / norm(subgradient, 2) ^ 2;
    
    
    currentPoint = currentPoint + alpha * subgradient;

    if strcmpi( options.verbose, 'iter')
        fprintf('I: %d, D: %f, P: %f, T: %f, stepSize: %f, nonZeroSubgr: %d\n', iIteration, dualPlot(iIteration),  primalPlot( iIteration ), timePlot(iIteration), alpha, sum( subgradient ~= 0));
    end
    
    % check stopping criteria
    
    if sum(subgradient .^ 2) == 0 
        reachedMaxIter = false;
        break;
    end
    
    if abs(bestPrimalValue - bestDualValue) < options.funcEps
        reachedMaxIter = false;
        break;
    end
    
    if  norm( alpha * subgradient ) < options.argTolEps
        reachedMaxIter = false;
        break;
    end
    
    if toc(tStart) - timeGarbage >= options.maxTime
        warning('maximizeSubgradientAdaptive:timeTermination','Process terminated: time elapsed\n');
        reachedMaxIter = false;
        break;        
    end
end

if reachedMaxIter
    warning('maximizeSubgradientAdaptive:maxIterTermination','Process terminated: maximum number of iterations reached\n');
end

lastNonNan = find(~isnan(timePlot), 1, 'last');
timePlot = timePlot(1 : lastNonNan);
dualPlot = dualPlot(1 : lastNonNan);
primalPlot = primalPlot(1 : lastNonNan);

end

function s = setDefaultField(s, f, value)
    if ~isfield(s, f)
        s.(f) = value;
    end
end

