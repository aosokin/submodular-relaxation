function [dualValue, subgradient, primalLabeling] = computeSmrDual_pairwisePotts(dataCost, neighbors, dualVars)
%computeSmrDual_pairwisePotts computes the value of the dual function in SMR method for pairwise energy with Potts potentials
%
% The function minimizes the Lagrangian over binary variables Y given duals variables D:
% L(Y, D) = \sum_i \sum_p U_{ip} y_{ip} + \sum_{ij} P_{ij} \sum_{p} 0.5 * ( [ y_{ip} == 1][y_{ip} == 0] + [ y_{ip} == 0][y_{ip} == 1] ) ...
%     +  \sum_i d_i ( \sum_p y_{ip} - 1)
%
% [dualValue, subgradient, primalLabeling]= computeSmrDual_pairwisePotts(dataCost, neighbors, dualVars)
%
% INPUT
%   dataCost   - unary potentials ( double[ numLabels x numNodes ])
%   neighbors  - paiwise Potts potentials ( sparse double[ numNodes x numNodes ]). 
%           The function uses only upper triangle of this matrix. All entries have to be non-negative.
%   dualVars   - vector of dual varuables ( double[ numNodes x 1 ])
%
% OUTPUT
%   dualValue - the value of the dual function
%   subgradient - value of subgradient
%   primalLabeling - the estimate of primal labeling
%
%   Depends on mexWrappers/graphCutMex_BoykovKolmogorov
%
% Anton Osokin (firstname.lastname@gmail.com),  24.09.2014

% check the input
if ~isnumeric(dataCost) || ~ismatrix(dataCost)
    error('computeSmrDual_pairwisePotts:badDataCost', 'dataCost should be a matrix  numLabels x numNodes');
end
dataCost = double(dataCost);
numNodes = size(dataCost, 2);
numLabels = size(dataCost, 1);

if ~isnumeric(neighbors) || ~ismatrix(neighbors) || ~issparse(neighbors) || size(neighbors, 1) ~= numNodes || size(neighbors, 2) ~= numNodes
    error('computeSmrDual_pairwisePotts:badNeighbors', 'neighbors should be a sparse matrix numNodes x numNodes');
end

if ~isnumeric(dualVars) || ~iscolumn(dualVars) || length(dualVars) ~= numNodes
    error('computeSmrDual_pairwisePotts:badDualVars', 'dualVars should be a column vector of length numNodes');
end
dualVars = double(dualVars);

subEnergy = nan(numLabels, 1);
labelsQp = nan(numNodes, numLabels);

% construct edges for a graph cut
[rowNeighbor, colNeighbor, weightNeighbor] = find(neighbors);
deleteMask = rowNeighbor >= colNeighbor;
rowNeighbor( deleteMask ) = [];
colNeighbor( deleteMask ) = [];
weightNeighbor( deleteMask ) = [];
nonTermEdgesWeights = [rowNeighbor, colNeighbor, 0.5 * weightNeighbor, 0.5 * weightNeighbor];

% construct unary terms for a graph cut
termEdgeWeight  = dataCost';

% run graph cuts
for iLabel = 1 : numLabels
    [subEnergy(iLabel), labelsQp(:, iLabel)] = graphCutMex([termEdgeWeight(:, iLabel) + dualVars, zeros(numNodes, 1)], nonTermEdgesWeights);
end
dualValue = sum(subEnergy) - sum(dualVars);

% get the primal estimate
if nargout > 2
    [~, primalLabeling] = max(labelsQp, [], 2);
end

%Compute subgradient
subgradient = sum(labelsQp, 2) - 1;

end
