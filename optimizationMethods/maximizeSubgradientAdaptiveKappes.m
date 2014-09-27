function [ bestPoint, bestVal, bestPrimal, bestLabeling, timePlot, funcPlot, primalPlot ] = maximizeSubgradientAdaptiveKappes( numVars, oracle, options )
% maximizeSubgradientAdaptive maximizes the dual function using adaptive subgradient method
%
% [ bestPoint, bestVal, timePlot, funcPlot, primalPlot ] = maximizeSubgradientAdaptive( oracle, options );
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
options = setDefaultField(options, 'epsScale', 1);
options = setDefaultField(options, 'initDelta', 1e+2);
options = setDefaultField(options, 'verbose', 'iter');
options = setDefaultField(options, 'argTolEps', 1e-20);
options = setDefaultField(options, 'gamma', 1.0);
options = setDefaultField(options, 'funcEps', 1e-5);



bestPoint = nan(numVars, 1);
bestVal = -inf;
bestPrimal = inf;
bestLabeling = [];
timePlot = nan(options.maxIter, 1);
funcPlot = nan(options.maxIter, 1);
primalPlot = nan(options.maxIter, 1);

tStart = tic;
timeGarbage = 0;
numOracleCalls = 0;
reachedMaxIter = true;

if isempty(options.initialPoint)
    initialPoint = zeros(numVars, 1);
else
    initialPoint = options.initialPoint;
end

currentPoint = initialPoint;
delta = options.initDelta;

for iIteration = 1 : options.maxIter
    
    [fValue, subgradient, primalEstimate] = oracle( currentPoint );
    funcPlot(iIteration) = fValue;
    timePlot(iIteration) = toc(tStart) - timeGarbage;
    numOracleCalls = numOracleCalls + 1;
    % exclude primal estimation from time 
    tTemp = tic;
    [primalPlot( iIteration ), newLabeling ] = options.funcGetPrimalLabeling(primalEstimate);
%     timeGarbage = timeGarbage + toc(tTemp);
    
    if fValue > bestVal
        bestVal = fValue;
        bestPoint = currentPoint;
    end
    if primalPlot( iIteration ) < bestPrimal
        bestPrimal = primalPlot( iIteration );
        bestLabeling = newLabeling;
    end
    
    
   alpha = options.gamma * (bestPrimal - fValue) / norm(subgradient, 2) ^ 2;
    
    
    currentPoint = currentPoint + alpha * subgradient;

    if strcmpi( options.verbose, 'iter')
        fprintf('I: %d, D: %f, P: %f, T: %f, stepSize: %f, nonZeroSubgr: %d\n', iIteration, funcPlot(iIteration),  primalPlot( iIteration ), timePlot(iIteration), alpha, sum( subgradient ~= 0));
    end
    
    % check stopping criteria
    
    if sum(subgradient .^ 2) == 0 
        reachedMaxIter = false;
        break;
    end
    
    if abs(bestPrimal - bestVal) < options.funcEps
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
    if numOracleCalls >= options.maxOracleCalls
        warning('maximizeSubgradientAdaptive:oracleTermination','Process terminated: maximum number of oracle calls reached\n');
        reachedMaxIter = false;
        break;        
    end
end

if reachedMaxIter
    warning('maximizeSubgradientAdaptive:maxIterTermination','Process terminated: maximum number of iterations reached\n');
end

end

function s = setDefaultField(s, f, value)
    if ~isfield(s, f)
        s.(f) = value;
    end
end

