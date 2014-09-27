function [ bestPoint, bestVal, bestPrimal, bestLabeling, timePlot, funcPlot, primalPlot ] = maximizeBundleMethodAggregated( numVars, oracle, options )
% maximizeBundleMethodAggregated maximizes the dual function using the bundle method with aggregation scheme by Kiwiel (1983)
%
% [ bestPoint, bestVal, timePlot, funcPlot, primalPlot ] = maximizeBundleMethodAggregated( oracle, options );
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
options = setDefaultField(options, 'nullSeriosStepRatio', 0.1); % m_r
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

[fValue, subgradient] = oracle( initialPoint );
funcPlot(1) = fValue;
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
    funcPlot(iIteration + 1) = fYk1;
    timePlot(iIteration + 1) = toc(tStart) - timeGarbage;
    numOracleCalls = numOracleCalls + 1;
    % exclude primal estimation from time 
    tTemp = tic;
    [primalPlot( iIteration + 1 ), newLabeling ] = options.funcGetPrimalLabeling(primalEstimate);
%     timeGarbage = timeGarbage + toc(tTemp);
    
    if fYk1 > bestVal
        bestVal = fYk1;
        bestPoint = y_k1;
    end
    if primalPlot( iIteration + 1 ) < bestPrimal
        bestPrimal = primalPlot( iIteration + 1 );
        bestLabeling = newLabeling;
    end
    
    % choose next stepSize
    nextInverseStepSize = 1 / (options.gammaStepSize * (min( bestPrimal, options.primalEstimate ) - bestVal) / norm(subgradient));
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
        fprintf('I: %d, S: %s, D: %f, P: %f, T: %f, invStepSize: %f, ratio: %f\n', iIteration, stepType, funcPlot(iIteration + 1), primalPlot( iIteration + 1), timePlot(iIteration + 1), inverseStepSize, (fYk1 - fXk) / v_k);
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







