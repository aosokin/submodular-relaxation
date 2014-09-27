function [funcValue, subgradient, primalLabeling] = computeNsmrDual_pairwisePotts(dataCost, neighbors, dualVars)
% computeNsmrDual_pairwisePotts computes the value of the dual function in SMD method for pairwise energy with Potts potentials
%
% The function minimizes the NSMR Lagrangian over binary variables Y given duals variables D:
% L(Y, D) = \sum_i \sum_p U_{ip} y_{ip} + \sum_{ij} P_{ij} \sum_{p} 0.5 * ( [ y_{ip} == 1][y_{ip} == 0] + [ y_{ip} == 0][y_{ip} == 1] ) ...
%     +  \sum_i d_i ( \sum_p y_{ip} - 1)
%
% [funcValue, subgradient, primalLabeling]= computeNsmrDual_pairwisePotts(dataCost, neighbors, dualVars)
%
% INPUT
%   dataCost   - unary potentials ( double[ numLabels x numNodes ])
%   neighbors  - paiwise Potts potentials ( sparse double[ numNodes x numNodes ]).
%       The function uses only upper triangle of this matrix. Can handle negative entries.
%   dualVars   - vector of dual varuables ( double[ numNodes x 1 ])
%
% OUTPUT
%   funcValue - the value of the dual function
%   subgradient - value of subgradient
%   primalLabeling - the estimate of primal labeling
%
% Anton Osokin (firstname.lastname@gmail.com),  24.09.2014

% check for the mex-file
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

subLowerBound = nan(numLabels, 1);
labelsQp = nan(numNodes, numLabels);

% construct edges for a graph cut
[rowNeighbor, colNeighbor, weightNeighbor] = find(neighbors);
deleteMask = rowNeighbor >= colNeighbor;
rowNeighbor( deleteMask ) = [];
colNeighbor( deleteMask ) = [];
weightNeighbor( deleteMask ) = [];
nEdges = length(weightNeighbor);
pairwiseTerms = [rowNeighbor, colNeighbor, zeros(nEdges, 1), 0.5 * weightNeighbor, 0.5 * weightNeighbor, zeros(nEdges, 1)];

% construct unary terms for a graph cut
termEdgeWeight  = dataCost';

for iLabel = 1 : numLabels
    unaryTerms = [zeros(numNodes, 1), termEdgeWeight(:, iLabel) + dualVars];
    [subLowerBound(iLabel), labelsQp(:, iLabel)] = qpboMex(unaryTerms, pairwiseTerms);
end
funcValue = sum(subLowerBound) - sum(dualVars);
labelsQp(labelsQp < 0) = 0.5;

% get the primal estimate
if nargout > 2
    [~, primalLabeling] = max(labelsQp, [], 2);
end

%Compute subgradient
subgradient = sum(labelsQp, 2) - 1;

end
