function [ bestPoint, bestVal, bestPrimal, bestLabeling, timePlot, funcPlot, primalPlot ] = maximizeBundleMethodFixedSize( numVars, oracle, options )
% maximizeBundleMethodFixedSize maximizes the dual function using the bundle method with the bundle of the fixed size
%
% [ bestPoint, bestVal, timePlot, funcPlot, primalPlot ] = maximizeBundleMethodFixedSize( oracle, options );
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
options = setDefaultField(options, 'maxBundleSize', 10);
options = setDefaultField(options, 'verbose', 'iter');
options = setDefaultField(options, 'nullSeriosStepRatio', 0.1); % m_L
options = setDefaultField(options, 'gammaStepSize', 0.1); % \gamma
options = setDefaultField(options, 'initialInverseStepSize', 1e+2); % w_0
options = setDefaultField(options, 'minInverseStepSize', 1e-10); % w_{min}
options = setDefaultField(options, 'maxInverseStepSize', 1e+7); % w_{max}
options = setDefaultField(options, 'primalEstimate', inf); 
options = setDefaultField(options, 'nullSeriosStepThreshhold', 5); 
options = setDefaultField(options, 'nullSeriosStepDivider', 2); 


bestPoint = nan(numVars, 1);
bestVal = -inf;
bestPrimal = inf;
bestLabeling = [];
timePlot = nan(options.maxIter, 1);
funcPlot = nan(options.maxIter, 1);
primalPlot = nan(options.maxIter, 1);

tStart = tic;
timeGarbage = 0;
numOracleCalls = 1;
reachedMaxIter = true;


if isempty(options.initialPoint)
    initialPoint = -(1e-3)*ones(numVars, 1);
else
    initialPoint = options.initialPoint;
end

currentPoint = initialPoint;
currentMeanPoint = currentPoint;
[fValue, subgradient] = oracle( currentPoint );
funcMeanValue = fValue;
funcPlot(1) = fValue;
timePlot(1) = toc(tStart); 

currentBundle = struct;
currentBundle.maxBundleSize = options.maxBundleSize;
currentBundle.numDimension = numVars;
currentBundle.size = 0;
% currentBundle.points = zeros(currentBundle.numDimension, currentBundle.maxBundleSize); % we do not actually need this
currentBundle.fValues = zeros(currentBundle.maxBundleSize, 1);
currentBundle.subgradients = single( zeros(currentBundle.numDimension, currentBundle.maxBundleSize) );
currentBundle.pointSubgradProd = zeros(currentBundle.maxBundleSize, 1);
currentBundle.grammMatrix = zeros(currentBundle.maxBundleSize, currentBundle.maxBundleSize);


inverseStepSize = options.initialInverseStepSize;
nullPositiveStepCounter = 0;

for iIteration = 1 : options.maxIter

    % make a step 
    currentBundle = updateBundle( currentBundle, currentPoint, fValue, subgradient );
    
    [nextPoint, bundleValue] = bundleProxStep( currentBundle, currentMeanPoint, inverseStepSize);
    
    % call the oracle
    [fValue, subgradient, primalEstimate] = oracle( nextPoint );
    funcPlot(iIteration + 1) = fValue;
    timePlot(iIteration + 1) = toc(tStart) - timeGarbage;
    numOracleCalls = numOracleCalls + 1;
    % exclude primal estimation from time 
    tTemp = tic;
    [primalPlot( iIteration + 1 ), newLabeling ]= options.funcGetPrimalLabeling(primalEstimate);
%     timeGarbage = timeGarbage + toc(tTemp);
    
    if fValue > bestVal
        bestVal = fValue;
        bestPoint = currentPoint;
    end
    if primalPlot( iIteration + 1 ) < bestPrimal
        bestPrimal = primalPlot( iIteration + 1 );
        bestLabeling = newLabeling;
    end
    
    % compute the next inv step size
    nextInverseStepSize = 1 / (options.gammaStepSize * (min( bestPrimal, options.primalEstimate )  - bestVal) / norm(subgradient));
    nextInverseStepSize = max(nextInverseStepSize, options.minInverseStepSize);
    nextInverseStepSize = min(nextInverseStepSize, options.maxInverseStepSize);
    
    % decide between a serious and a null step
    seriousStepRatio = (fValue - funcMeanValue) / (bundleValue - funcMeanValue);
    if (fValue > funcMeanValue && seriousStepRatio > options.nullSeriosStepRatio)
        % serious step
        currentMeanPoint = nextPoint;
        currentPoint = nextPoint;
        funcMeanValue = fValue;
        inverseStepSize = nextInverseStepSize;
        stepType = 'SERIOS';
        nullPositiveStepCounter = 0;
    else 
        % null step
        currentPoint = nextPoint;
        % Bogdan says not to change inverseStepSize at a null step
%         inverseStepSize = nextInverseStepSize;
        stepType = 'NULL';
        if seriousStepRatio > 0
            nullPositiveStepCounter = nullPositiveStepCounter + 1;
        end
    end
    
    if nullPositiveStepCounter > options.nullSeriosStepThreshhold
        nullPositiveStepCounter = 0;
        options.nullSeriosStepRatio = options.nullSeriosStepRatio / options.nullSeriosStepDivider;
    end

    if strcmpi( options.verbose, 'iter')
        fprintf('I: %d, S: %s, D: %f, T: %f, invStepSize: %f, ratio: %f\n', iIteration, stepType, funcPlot(iIteration + 1), timePlot(iIteration + 1), inverseStepSize, seriousStepRatio);
    end
    
    % move to next iteration
%      inverseStepSize = options.initialInverseStepSize;
    
    % check stopping criteria

    if sum( subgradient.^2 ) == 0
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







