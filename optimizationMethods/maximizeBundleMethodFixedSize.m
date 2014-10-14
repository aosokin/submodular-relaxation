function [ bestDualPoint, bestDualValue, bestPrimalValue, bestPrimalLabeling, timePlot, dualPlot, primalPlot ] = maximizeBundleMethodFixedSize( numVars, oracle, options )
%maximizeBundleMethodFixedSize maximizes the dual function using the bundle method 
% with bundle of limited size.
%
% Please cite the following paper if you use this code:
%   J. H. Kappes, B. Savchynskyy, and C. Schnoerr, 
%   A bundle approach to efficient MAP-inference by Lagrangian relaxation,
%   CVPR, 2012.
%
% Usage:
% [ bestDualPoint, bestDualValue, bestPrimalValue, bestPrimalLabeling, timePlot, dualPlot, primalPlot ] = maximizeBundleMethodFixedSize( numVars, oracle, options );
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
%       maxBundleSize - the maximum size of the bundle (default: 10);
%       nullSeriosStepRatio - the method parameter m_L; see the journal paper for the details; (default: 0.1)
%       gammaStepSize - the method parameter \gamma; see the journal paper for the details; (default: 0.1)
%       initialInverseStepSize - the method parameter w_0; see the journal paper for the details; (default: 1e+2)
%       minInverseStepSize - the method parameter w_{min}; see the journal paper for the details; (default: 1e-10)
%       maxInverseStepSize - the method parameter w_{max}; see the journal paper for the details; (default: 1e+7)
%       primalEstimate - the initial estimate of the primal function; (default: inf)
%       nullSeriosStepThreshhold - after nullSeriosStepThreshhold number of null-steps in a row the nullSeriosStepRatio is divided by nullSeriosStepDivider; (default: 5)
%       nullSeriosStepDivider - after nullSeriosStepThreshhold number of null-steps in a row the nullSeriosStepRatio is divided by nullSeriosStepDivider; (default: 2)
%
% OUTPUT:
%       bestDualPoint - the value of the variables that delivers the maximum-found function value
%       bestDualValue - the function value at bestDualPoint
%       bestPrimalValue - the lowest-found energy
%       bestPrimalLabeling - the labeling that delivers bestPrimalValue
%       timePlot, dualPlot, primalPlot - running time, the dual and the primal values for each oracle call 
%
%   Depends on optimizationMethods/bundleMethod
%
% Anton Osokin (firstname.lastname@gmail.com),  16.05.2013

if ~exist( 'options', 'var');
    options = struct;
end
options = setDefaultField(options, 'maxIter', 100);
options = setDefaultField(options, 'maxTime', inf);
options = setDefaultField(options, 'funcGetPrimalLabeling', @(x) funcGetPrimalLabeling_default(x));
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


bestDualPoint = nan(numVars, 1);
bestDualValue = -inf;
bestPrimalValue = inf;
bestPrimalLabeling = [];
timePlot = nan(options.maxIter, 1);
dualPlot = nan(options.maxIter, 1);
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
dualPlot(1) = fValue;
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
    dualPlot(iIteration + 1) = fValue;
    timePlot(iIteration + 1) = toc(tStart) - timeGarbage;
    numOracleCalls = numOracleCalls + 1;
    % exclude primal estimation from time 
    tTemp = tic;
    [primalPlot( iIteration + 1 ), newLabeling ]= options.funcGetPrimalLabeling(primalEstimate);
%     timeGarbage = timeGarbage + toc(tTemp);
    
    if fValue > bestDualValue
        bestDualValue = fValue;
        bestDualPoint = currentPoint;
    end
    if primalPlot( iIteration + 1 ) < bestPrimalValue
        bestPrimalValue = primalPlot( iIteration + 1 );
        bestPrimalLabeling = newLabeling;
    end
    
    % compute the next inv step size
    nextInverseStepSize = 1 / (options.gammaStepSize * (min( bestPrimalValue, options.primalEstimate )  - bestDualValue) / norm(subgradient));
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
        fprintf('I: %d, S: %s, D: %f, T: %f, invStepSize: %f, ratio: %f\n', iIteration, stepType, dualPlot(iIteration + 1), timePlot(iIteration + 1), inverseStepSize, seriousStepRatio);
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







