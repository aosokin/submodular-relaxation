function [ bestDualPoint, bestDualValue, bestPrimalValue, bestPrimalLabeling, timePlot, dualPlot, primalPlot ] = maximizeBundleMethodAggregated( numVars, oracle, options )
%maximizeBundleMethodAggregated maximizes the dual function using the bundle method 
% with bundle of limited size.
%
% Please cite the following papers if you use this code:
%   1) J. H. Kappes, B. Savchynskyy, and C. Schnoerr, 
%   A bundle approach to efficient MAP-inference by Lagrangian relaxation,
%   CVPR, 2012.
%   2) K. Kiwiel, An aggregate subgradient method for nonsmooth convex minimization,
%   Mathematical Programming, vol. 27, pp. 320–341, 1983.
%
% Usage:
% [ bestDualPoint, bestDualValue, bestPrimalValue, bestPrimalLabeling, timePlot, dualPlot, primalPlot ] = maximizeBundleMethodAggregated( numVars, oracle, options );
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
%       verbose - the verbosity level, could be 'iter', or 'none' (default: 'iter')
%       nullSeriosStepRatio - the method parameter m_r; see the journal paper for the details; (default: 0.1)
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
% Anton Osokin (firstname.lastname@gmail.com),  16.05.2013

if ~exist( 'options', 'var');
    options = struct;
end
options = setDefaultField(options, 'maxIter', 100);
options = setDefaultField(options, 'maxTime', inf);
options = setDefaultField(options, 'funcGetPrimalLabeling', @(x) funcGetPrimalLabeling_default(x));
options = setDefaultField(options, 'initialPoint', []);
options = setDefaultField(options, 'verbose', 'iter');
options = setDefaultField(options, 'nullSeriosStepRatio', 0.1); % m_r
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

[fValue, subgradient] = oracle( initialPoint );
dualPlot(1) = fValue;
timePlot(1) = toc(tStart);

aggrBundle = subgradient;
alpha_k = 0;
eps_k = 0;
nu_k = nan;
fYk1 = nan;
fXk = fValue;
x_k = initialPoint;
y_k = nan;
x_k1 = nan;
y_k1 = nan;
nullStepCounter = 0;
inverseStepSize = options.initialInverseStepSize;


for iIteration = 1 : options.maxIter

    % make a step 

    if  aggrBundle' * aggrBundle + subgradient' * subgradient - 2 * aggrBundle' * subgradient ~= 0
        mu = ( eps_k - alpha_k - aggrBundle' * subgradient + aggrBundle' * aggrBundle ) / ( aggrBundle' * aggrBundle + subgradient' * subgradient - 2 * aggrBundle' * subgradient);
    else
        mu = nan;
    end
    
    func = @(curMu) (0.5 * norm( (1 - curMu) * aggrBundle + curMu * subgradient, 2 )^2 + (1 - curMu) * eps_k + curMu * alpha_k);
    
    curVal = func(mu);
    if func(0) < curVal || isnan(mu)
        curVal = func(0);
        mu = 0;
    end
    if func(1) < curVal
        curVal = func(1);
        mu = 1;
    end

    
    aggrBundle = (1 - mu) * aggrBundle + mu * subgradient;
    nu_k = (1 - mu) * eps_k + mu * alpha_k;
    v_k = ( norm( aggrBundle, 2 )^2 + nu_k );
    
    %t_Rk = 1e-3;
    t_Rk = 1 / inverseStepSize;
    y_k1 = x_k + aggrBundle * t_Rk;
    
    % call the oracle
    [fYk1, subgradient, primalEstimate] = oracle( y_k1 );
    dualPlot(iIteration + 1) = fYk1;
    timePlot(iIteration + 1) = toc(tStart) - timeGarbage;
    numOracleCalls = numOracleCalls + 1;
    % exclude primal estimation from time 
    tTemp = tic;
    [primalPlot( iIteration + 1 ), newLabeling ] = options.funcGetPrimalLabeling(primalEstimate);
%     timeGarbage = timeGarbage + toc(tTemp);
    
    if fYk1 > bestDualValue
        bestDualValue = fYk1;
        bestDualPoint = y_k1;
    end
    if primalPlot( iIteration + 1 ) < bestPrimalValue
        bestPrimalValue = primalPlot( iIteration + 1 );
        bestPrimalLabeling = newLabeling;
    end
    
    % choose next stepSize
    nextInverseStepSize = 1 / (options.gammaStepSize * (min( bestPrimalValue, options.primalEstimate ) - bestDualValue) / norm(subgradient));
    nextInverseStepSize = max(nextInverseStepSize, options.minInverseStepSize);
    nextInverseStepSize = min(nextInverseStepSize, options.maxInverseStepSize);
    
    % choose null or serious step
    if fYk1 >= fXk + options.nullSeriosStepRatio * v_k
        % serious step
        % t_Lk = t_Rk;
        x_k1 = y_k1;
        fXk1 = fYk1;
        stepType = 'SERIOUS';
        inverseStepSize = nextInverseStepSize;
        nullStepCounter = 0;
    else
        % null step 
        % t_Lk = 0;
        x_k1 = x_k;
        fXk1 = fXk;
        stepType = 'NULL';
        if fYk1 >= fXk
            nullStepCounter = nullStepCounter + 1;
        end
    end
    alpha_k = -( fXk1 - fYk1 - subgradient' * ( x_k1 - y_k1 ) );
    eps_k = nu_k -( fXk1 - fXk - aggrBundle' * ( x_k1 - x_k ) );

    if (nullStepCounter > options.nullSeriosStepThreshhold)
        options.nullSeriosStepRatio = options.nullSeriosStepRatio / options.nullSeriosStepDivider;
    end
        
    if strcmpi( options.verbose, 'iter')
        fprintf('I: %d, S: %s, D: %f, P: %f, T: %f, invStepSize: %f, ratio: %f\n', iIteration, stepType, dualPlot(iIteration + 1), primalPlot( iIteration + 1), timePlot(iIteration + 1), inverseStepSize, (fYk1 - fXk) / v_k);
    end
    
    % move to next iteration
    x_k = x_k1;
    x_k1 = nan;
    fXk = fXk1;
    fXk1 = nan;
    y_k1 = nan;
    fYk1 = nan;

    
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







