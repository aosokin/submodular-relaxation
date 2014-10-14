function [dualValue, subgradient, primalLabeling] = computeSmrDualSubtraction_pairwisePotts(dataCost, neighbors, dualVars)
% computeSmrDualSubtraction_pairwisePotts computes the value of the dual function in SMD method for pairwise energy with Potts potentials
%
% The function minimizes the SMR Lagrangian (subtraction trick) over binary variables Y given duals variables D:
% L(Y, D) = \sum_i \sum_p U_{ip} y_{ip} + \sum_{ij} 0.5 * P_{ij} \sum_{p \neq q} ( [ y_{ip} == 1][y_{iq} == 0] + [ y_{ip} == 0][y_{iq} == 1] ) ...
%     +  \sum_i d_i ( \sum_p y_{ip} - 1)
%
% [dualValue, subgradient, primalLabeling]= computeSmrDualSubtraction_pairwisePotts(dataCost, neighbors, dualVars)
%
% INPUT
%   dataCost   - unary potentials ( double[ numLabels x numNodes ])
%   neighbors  - paiwise Potts potentials ( sparse double[ numNodes x numNodes ]).
%       The function uses only upper triangle of this matrix. Can handle negative entries.
%   dualVars   - vector of dual varuables ( double[ numNodes x 1 ])
%
% OUTPUT
%   dualValue - the value of the dual function
%   subgradient - value of subgradient
%   primalLabeling - the estimate of primal labeling
%
% Depends on mexWrappers/graphCutMex_BoykovKolmogorov
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

% construct unary terms for a graph cut
dataCostLambda = bsxfun(@plus, dataCost, dualVars');

% construct edges for a graph cut
[rowNeighbor, colNeighbor, weightNeighbor] = find(neighbors);
deleteMask = rowNeighbor >= colNeighbor;
rowNeighbor( deleteMask ) = [];
colNeighbor( deleteMask ) = [];
weightNeighbor( deleteMask ) = [];

positiveEdges = find(weightNeighbor > 0);
negativeEdges = find(weightNeighbor < 0);

pairwiseTerms = nan( length(negativeEdges) * numLabels * (numLabels - 1) + ...
                     length(positiveEdges) * numLabels, 4);

nEdgesSet = 0;

% add positive edges
for iEdge = 1 : length(positiveEdges)
    curRow = rowNeighbor( positiveEdges(iEdge) );
    curCol = colNeighbor( positiveEdges(iEdge) );
    curWeight = weightNeighbor( positiveEdges(iEdge) );
    
    newWeights = 0.5 * curWeight * ones(numLabels, 1);
    id1 = (1 : numLabels)' + (curRow - 1) * numLabels;
    id2 = (1 : numLabels)' + (curCol - 1) * numLabels;
    pairwiseTerms(nEdgesSet + 1 : nEdgesSet + numLabels, :) = ...
        [ id1, id2, newWeights, newWeights ];
    
    nEdgesSet = nEdgesSet + numLabels;
end

% add negative edges
id1 = repmat(1 : numLabels, [numLabels, 1]);
id2 = repmat((1 : numLabels)', [1, numLabels]);
idMask = id1 ~= id2;
id1_b = id1(idMask(:));
id2_b = id2(idMask(:));
for iEdge = 1 : length(negativeEdges)
    curRow = rowNeighbor( negativeEdges(iEdge) );
    curCol = colNeighbor( negativeEdges(iEdge) );
    curWeight = weightNeighbor( negativeEdges(iEdge) );
    
    dataCostLambda(:, curRow) = dataCostLambda(:, curRow) + curWeight * 0.5 * (numLabels - 1);
    dataCostLambda(:, curCol) = dataCostLambda(:, curCol) + curWeight * 0.5 * (numLabels - 1);
    
    newWeights = -0.5 * curWeight * ones( numLabels * (numLabels - 1), 1);
    
    id1 = id1_b + (curCol - 1) * numLabels;
    id2 = id2_b + (curRow - 1) * numLabels;

%     % slow
%     id1 = bsxfun(@plus, 1 : numLabels, (curCol - 1) * numLabels * ones(numLabels, 1));
%     id2 = bsxfun(@plus, (1 : numLabels)', (curRow - 1) * numLabels * ones(1, numLabels));
%     id1 = id1(idMask);
%     id2 = id2(idMask);
    
    pairwiseTerms(nEdgesSet + 1 : nEdgesSet + numLabels * (numLabels - 1), :) = ...
        [ id1(:), id2(:), newWeights, newWeights ];
    
    nEdgesSet = nEdgesSet + numLabels * (numLabels - 1);
end

termWeights  = [dataCostLambda(:), zeros(numNodes * numLabels, 1)];

[cutValue, cutLabels] = graphCutMex(termWeights, pairwiseTerms);
labelsQp = reshape( cutLabels, [numLabels, numNodes] )';

dualValue = cutValue - sum(dualVars); 
% get the primal estimate
if nargout > 2
    [~, primalLabeling] = max(labelsQp, [], 2);
end

global globalNumOracleCalls

%Compute subgradient
subgradient = sum(labelsQp, 2) - 1;

i = 1;
    
end
