function [dualValue, subgradient, primalLabeling] = computeSmrDualDynamic_pairwisePotts(dataCost, neighbors, dualVars)
% computeSmrDualDynamic_pairwisePotts computes the value of the SMR dual function in Smr method for pairwise energy with Potts potentials
%
% The function minimizes the Lagrangian w.r.t. binary variables  Y given duals variables D:
% L(Y, D) = \sum_i \sum_p U_{ip} y_{ip} + \sum_{ij} P_{ij} \sum_{p} 0.5 * ( [ y_{ip} == 1][y_{ip} == 0] + [ y_{ip} == 0][y_{ip} == 1] ) ...
%     +  \sum_i d_i ( \sum_p y_{ip} - 1)
%
%   This function makes use of dynamic graph cuts to compute the updates faster
%   the following global variables are used: computeSmrDualDynamic_pairwisePotts_graphHandle, computeSmrDualDynamic_pairwisePotts_lastPoint
%
% [dualValue, subgradient, primalLabeling]= computeSmrDualDynamic_pairwisePotts(dataCost, neighbors, dualVars)
%
% INPUT
%   dataCost   - unary potentials ( double[ numLabels x numNodes ])
%   neighbors  - paiwise Potts potentials ( sparse double[ numNodes x numNodes ]). The function uses only upper triangle of this matrix.
%   dualVars   - vector of dual varuables ( double[ numNodes x 1 ])
%
% OUTPUT
%   dualValue - the value of the dual function
%   subgradient - value of subgradient
%   primalLabeling - the estimate of primal labeling
%
% CAUTION! do not forget to call computeSmrDualDynamic_pairwisePotts_clearGlobal after the optimization is finished
%
% Depends on mexWrappers/graphCutDynamicMex_BoykovKolmogorov
%
% Anton Osokin (firstname.lastname@gmail.com),  24.09.2014

% check the input
if ~isnumeric(dataCost) || ~ismatrix(dataCost)
    error('computeSmrDualDynamic_pairwisePotts:badDataCost', 'dataCost should be a matrix  numLabels x numNodes');
end
dataCost = double(dataCost);
numNodes = size(dataCost, 2);
numLabels = size(dataCost, 1);

if ~isnumeric(neighbors) || ~ismatrix(neighbors) || ~issparse(neighbors) || size(neighbors, 1) ~= numNodes || size(neighbors, 2) ~= numNodes
    error('computeSmrDualDynamic_pairwisePotts:badNeighbors', 'neighbors should be a sparse matrix numNodes x numNodes');
end

if ~isnumeric(dualVars) || ~iscolumn(dualVars) || length(dualVars) ~= numNodes
    error('computeSmrDualDynamic_pairwisePotts:badDualVars', 'dualVars should be a column vector of length numNodes');
end
dualVars = double(dualVars);

subEnergy = nan(numLabels, 1);
labelsQp = nan(numNodes, numLabels);

% after this number of runs recompute graph cut from scratch
dynamicCutRebuildNumber = 20;

global computeSmrDualDynamic_pairwisePotts_graphHandle
global computeSmrDualDynamic_pairwisePotts_lastPoint
global computeSmrDualDynamic_pairwisePotts_dynamicNumber 

if isempty(computeSmrDualDynamic_pairwisePotts_graphHandle) || isempty(computeSmrDualDynamic_pairwisePotts_lastPoint) || isempty(computeSmrDualDynamic_pairwisePotts_dynamicNumber)...
        || ~iscell(computeSmrDualDynamic_pairwisePotts_graphHandle) || numel( computeSmrDualDynamic_pairwisePotts_graphHandle ) ~= numLabels ...
        || ~iscolumn(computeSmrDualDynamic_pairwisePotts_lastPoint) || length( computeSmrDualDynamic_pairwisePotts_lastPoint ) ~=  numNodes ...
        || ~isscalar(computeSmrDualDynamic_pairwisePotts_dynamicNumber) || ~isnumeric(computeSmrDualDynamic_pairwisePotts_dynamicNumber) ...
        || mod( computeSmrDualDynamic_pairwisePotts_dynamicNumber, dynamicCutRebuildNumber) == 0
    % remove all graphs if left
    if ~isempty(computeSmrDualDynamic_pairwisePotts_graphHandle) 
        for iLabel = 1 : length(computeSmrDualDynamic_pairwisePotts_graphHandle);
            deleteGraphCutDynamicMex( computeSmrDualDynamic_pairwisePotts_graphHandle{iLabel} )
        end
    end
    
    % construct edges for a graph cut
    [rowNeighbor, colNeighbor, weightNeighbor] = find(neighbors);
    deleteMask = rowNeighbor >= colNeighbor;
    rowNeighbor( deleteMask ) = [];
    colNeighbor( deleteMask ) = [];
    weightNeighbor( deleteMask ) = [];
    nonTermEdgesWeights = [rowNeighbor, colNeighbor, 0.5 * weightNeighbor, 0.5 * weightNeighbor];

    % construct unary terms for a graph cut
    termEdgeWeight  = dataCost';
    
    % store a point
    computeSmrDualDynamic_pairwisePotts_lastPoint = dualVars;
    computeSmrDualDynamic_pairwisePotts_graphHandle = cell(numLabels, 1);
    computeSmrDualDynamic_pairwisePotts_dynamicNumber = 1;

    for iLabel = 1 : numLabels
        [subEnergy(iLabel), labelsQp(:, iLabel),...
            computeSmrDualDynamic_pairwisePotts_graphHandle{ iLabel }] = ...
            graphCutDynamicMex([termEdgeWeight(:, iLabel) + dualVars, zeros(numNodes, 1)], nonTermEdgesWeights);
    end
else
    pointDifference = dualVars - computeSmrDualDynamic_pairwisePotts_lastPoint;
    
    changeMask = pointDifference ~= 0;
    numChanges = sum( changeMask );
    
%     fprintf('Updated %f%% nodes \n', numChanges / numNodes * 100);
  
    unaryUpdate = [find(pointDifference), pointDifference( changeMask ), zeros( numChanges, 1 )];
    for iLabel = 1 : numLabels
        [subEnergy(iLabel), labelsQp(:, iLabel)] = updateUnaryGraphCutDynamicMex( computeSmrDualDynamic_pairwisePotts_graphHandle{ iLabel }, unaryUpdate );
    end
    computeSmrDualDynamic_pairwisePotts_lastPoint = dualVars;
    computeSmrDualDynamic_pairwisePotts_dynamicNumber = computeSmrDualDynamic_pairwisePotts_dynamicNumber + 1;
end

% compute the lower bound
dualValue = sum(subEnergy) - sum(dualVars);

% get the primal estimate
if nargout > 2
    [~, primalLabeling] = max(labelsQp, [], 2);
end

%Compute subgradient
subgradient = sum(labelsQp, 2) - 1;

% % check
% [old_funcValue, old_subgradient] = computeSmrDual_pairwisePotts(dataCost, neighbors, dualVars);
% 
% if abs( old_funcValue - funcValue ) > 1e-5
%     f = 5;
% end
% 
% if any( old_subgradient ~= subgradient)
%     fprintf('Difference in subgradient: %d \n', sum(old_subgradient ~= subgradient) );
% end

end
