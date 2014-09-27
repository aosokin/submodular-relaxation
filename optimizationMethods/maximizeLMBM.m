function [ bestPoint, bestVal, bestPrimal, bestLabeling, timePlot, funcPlot, primalPlot ] = maximizeLMBM( numVars, oracle, options )
% maximizeLMBM maximizes the dual function using LMBM method
%
% [ bestPoint, bestVal, timePlot, funcPlot, primalPlot ] = maximizeLMBM( oracle, options );
%
% INPUT:
%
%
% OUTPUT:
%
%
% Anton Osokin (firstname.lastname@gmail.com),  25.05.2013



if ~exist( 'options', 'var');
    options = struct;
end
options = setDefaultField(options, 'maxIter', 100);
options = setDefaultField(options, 'maxOracleCalls', 1e+8);
options = setDefaultField(options, 'maxTime', 1e+8);
options = setDefaultField(options, 'funcGetPrimalLabeling', @(x) nan);
options = setDefaultField(options, 'initialPoint', []);
options = setDefaultField(options, 'maxBundleSize', 10);
options = setDefaultField(options, 'verbose', 'iter'); %final or none
options = setDefaultField(options, 'method', 'lmbm'); % lmbm or lbfgs
options = setDefaultField(options, 'scale', 1.0); % scale for dual variables


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

global globalOraclePointer
globalOraclePointer = oracle;

global globalPrimalFuncPointer
globalPrimalFuncPointer = options.funcGetPrimalLabeling;

global globalIterPrintFlag
globalIterPrintFlag = false;

global globalDualVarsScale
globalDualVarsScale = options.scale;

if isempty(options.initialPoint)
    initialPoint = -(1e-3)*ones(numVars, 1);
else
    initialPoint = options.initialPoint;
end

if strcmpi(options.verbose, 'iter')
    lmbmPrtLevel = 1;
    globalIterPrintFlag = true;
elseif strcmpi(options.verbose, 'final')
    lmbmPrtLevel = 1;
elseif strcmpi(options.verbose, 'none')
    lmbmPrtLevel = 0;
else
    warning('maximizeLMBM:wrongOptionsVerbose', 'Unknown value of options.verbose, expected iter, final, or none');
end

lmbmMaxTime = options.maxTime;
lmbmNA = options.maxBundleSize;
lmbmMCU = [];
lmbmMC = [];
lmbmPPAR = nan(8, 1);
lmbmIPAR = nan(7, 1);
lmbmIPAR(2) = options.maxIter;
lmbmIPAR(3) = options.maxOracleCalls;
if strcmpi( options.method, 'lmbm')
    lmbmIPAR(6) = 0; % 0 - LMBM, 1 - LBFGS
elseif strcmpi( options.method, 'lbfgs')
    lmbmIPAR(1) = 0; % 0 - LMBM, 1 - LBFGS
else
    warning('maximizeLMBM:wrongOptionsMethod', 'Unknown value of options.method, expected lmbm of lbfgs none');
end
    
    


[lmbmPoint, lmbmVal, ~, ~, lmbmExitFlag] = lmbm_driver('maximizeLMBM_oracleFunction', ...
   initialPoint * globalDualVarsScale, numVars, lmbmPrtLevel, lmbmMaxTime, lmbmNA, lmbmMCU, lmbmMC, lmbmPPAR, lmbmIPAR); 
%[X,FVAL,NITER,NFEVAL,TERM,TIME]=LMBM_DRIVER(F,X0,N,PRINT,[MAXTIME,NA,MCU,MC,RPAR,IPAR]) 

fprintf('LMBM terminated with exit flag %d\n', lmbmExitFlag);

bestPoint = lmbmPoint;
bestVal = -lmbmVal;
if abs( globalBestVal - bestVal ) > 1e-5 * (max(abs(globalBestVal), abs(bestVal)) + 1e-5)
    warning('maximizeLMBM:maxLowerBoundMismatch', ['Best LBs computed by LMBM and explicitly do not match, diff: ', num2str( globalBestVal - bestVal )]);
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

clear global globalOraclePointer
clear global globalPrimalFuncPointer
clear global globalIterPrintFlag
clear global globalDualVarsScale
end

function s = setDefaultField(s, f, value)
    if ~isfield(s, f)
        s.(f) = value;
    end
end

% function [fValNegative, subgradientNegative] = inverseOracle(x, oracle, getPrimalLabeling)
%     [fVal, subgradient, primalEstimate] = oracle(x);
%     fValNegative = -fVal;
%     subgradientNegative = -subgradient;
%     
%     curStart = tic;
%     
%     global globalBestPoint
%     global globalBestVal
%     global globalBestPrimal
%     global globalBestLabeling
% 
%     global globalTimePlot
%     global globalFuncPlot
%     global globalPrimalPlot
% 
%     global globalTimeStart
%     global globalTimeGarbage
%     global globalNumOracleCalls
% 
%     globalNumOracleCalls = globalNumOracleCalls + 1;
%     globalFuncPlot( globalNumOracleCalls ) = fVal;
%     if fVal > globalBestVal
%         globalBestVal = fVal;
%         globalBestPoint = x;
%     end
%     
%     [curPrimal, newLabeling] = getPrimalLabeling( primalEstimate );
%     if curPrimal < globalBestPrimal
%         globalBestPrimal = curPrimal;
%         globalBestLabeling = newLabeling;
%     end
%     globalPrimalPlot( globalNumOracleCalls ) = curPrimal;
%     
%     globalTimePlot( globalNumOracleCalls ) = nan;
%     
%     curGarbage = toc( curStart );
%     globalTimeGarbage = globalTimeGarbage + curGarbage;
%     
%     globalTimePlot( globalNumOracleCalls ) = toc( globalTimeStart ) - globalTimeGarbage;
% end





