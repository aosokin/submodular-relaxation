function [ bestDualPoint, bestDualValue, bestPrimalValue, bestPrimalLabeling, timePlot, dualPlot, primalPlot ] = maximizeLMBM( numVars, oracle, options )
% maximizeLMBM maximizes the dual function with the LMBM method: http://napsu.karmitsa.fi/lmbm/
%
% Please cite the following paper if you use this code:
%   N. Haarala, K. Miettinen, and M. M. Makela, 
%   Globally convergent limited memory bundle method for large-scale nonsmooth optimization,
%   Mathematical Programming, vol. 109, no. 1, pp. 181–205, 2007.
%
% Usage:
% [ bestDualPoint, bestDualValue, bestPrimalValue, bestPrimalLabeling, timePlot, dualPlot, primalPlot ] = maximizeLMBM( numVars, oracle, options );
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
%       funcGetPrimalLabeling - function to compute the current estimate of the primal solution;
%           should be of format [ energy, labeling ] = funcGetPrimalLabeling( primalEstimate ); where
%               energy - is the energy value;
%               labeling - the current labeling;
%               primalEstimate - the initialization for the labeling coming from the dual
%           (default: nan for all labelings), consider using icmPottsMex
%       initialPoint - the initialization of the dual variables, (default: -(1e-3)*ones(numVars, 1)).
%       verbose - the verbosity level, could be 'iter', 'final', or 'none' (default: 'iter')
%       maxBundleSize - the maximum size of the bundle (default: 10);
%       method - 'lmbm' or 'lbfgs' (default: 'lmbm')
%       scale: scale of the dual variables (default: 1.0)
%
% OUTPUT:
%       bestDualPoint - the value of the variables that delivers the maximum-found function value
%       bestDualValue - the function value at bestDualPoint
%       bestPrimalValue - the lowest-found energy
%       bestPrimalLabeling - the labeling that delivers bestPrimalValue
%       timePlot, dualPlot, primalPlot - running time, the dual and the primal values for each oracle call 
%
%   Depends on LMBM: optimizationMethods/LMBM
%
% Anton Osokin (firstname.lastname@gmail.com),  16.05.2013


if ~exist( 'options', 'var');
    options = struct;
end
options = setDefaultField(options, 'maxIter', 100);
options = setDefaultField(options, 'funcGetPrimalLabeling', @(x) funcGetPrimalLabeling_default(x));
options = setDefaultField(options, 'initialPoint', []);
options = setDefaultField(options, 'maxBundleSize', 10);
options = setDefaultField(options, 'verbose', 'iter'); % 'iter', 'final' or 'none'
options = setDefaultField(options, 'method', 'lmbm'); % 'lmbm' or 'lbfgs'
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

lmbmMaxTime = 1e+8; %options.maxTime;
lmbmNA = options.maxBundleSize;
lmbmMCU = [];
lmbmMC = [];
lmbmPPAR = nan(8, 1);
lmbmIPAR = nan(7, 1);
lmbmIPAR(2) = options.maxIter;
lmbmIPAR(3) = 1e+8; %options.maxOracleCalls;
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

if lmbmPrtLevel ~= 0
    fprintf('LMBM terminated with exit flag %d\n', lmbmExitFlag);
end

bestDualPoint = lmbmPoint;
bestDualValue = -lmbmVal;
if abs( globalBestVal - bestDualValue ) > 1e-5 * (max(abs(globalBestVal), abs(bestDualValue)) + 1e-5)
    warning('maximizeLMBM:maxLowerBoundMismatch', ['Best LBs computed by LMBM and explicitly do not match, diff: ', num2str( globalBestVal - bestDualValue )]);
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

