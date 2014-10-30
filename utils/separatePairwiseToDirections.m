function [pairwiseDirectionalCosts, newNeighbors] = separatePairwiseToDirections( neighbors, gridSize )
% separatePairwiseToDirections extracts cost for vertical, horizontal, and diagonal edges from a sparse matrix
%
%   separatePairwiseToDirections is similar to separateVertHorCosts but with different interface
%
% 	Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

if ~isnumeric(gridSize) || ~iscolumn(gridSize) || length(gridSize) ~= 2 
    error('separatePairwiseToDirections:badGridSize', 'gridSize should be a column vector of length 2');
end
gridSize = double(gridSize);
numNodes = prod( gridSize );

if ~isnumeric(neighbors) || ~ismatrix(neighbors) || ~issparse(neighbors) || size(neighbors, 1) ~= numNodes || size(neighbors, 2) ~= numNodes
    error('separatePairwiseToDirections:badNeighbors', 'neighbors should be a sparse matrix numNodes x numNodes');
end


[neighRow, neighCol, neighWeight] = find(neighbors);
deleteMask = neighRow > neighCol;
neighRow( deleteMask ) = [];
neighCol( deleteMask ) = [];
neighWeight( deleteMask ) = [];

% [~, ids] = sortrows( [ neighRow, neighCol ] );
% neighRow = neighRow(ids);
% neighCol = neighCol(ids);
% neighWeight = neighWeight(ids);

vertMask = (neighRow + 1 == neighCol);
horMask = (neighRow + gridSize(1) == neighCol);
mainDiagMask = (neighRow + 1 + gridSize(1) == neighCol);
secondDiagMask = (neighRow - 1 + gridSize(1) == neighCol);

if sum( vertMask ) + sum( horMask ) + sum( mainDiagMask ) + sum( secondDiagMask ) ~= length(neighRow)
    warning('separateVertHorCosts:wrongNeighborhoodSystem', 'The graph is not 4-, 6-, or 8-connected grid!!!');
end

vertCost = reshape(neighWeight(vertMask), [gridSize(1) - 1, gridSize(2)]);
horCost = reshape(neighWeight(horMask), [gridSize(1), gridSize(2) - 1]);
if any(mainDiagMask(:))
    mainDiagCost = reshape(neighWeight(mainDiagMask), [gridSize(1) - 1, gridSize(2) - 1]);
end
if any(secondDiagMask(:))
    secondDiagCost = reshape(neighWeight(secondDiagMask), [gridSize(1) - 1, gridSize(2) - 1]);
end


toKeepMask = vertMask | horMask | mainDiagMask | secondDiagMask;
newNeighbors = sparse(neighRow(toKeepMask), neighCol(toKeepMask), neighWeight(toKeepMask), numNodes, numNodes);

pairwiseDirectionalCosts = struct;
pairwiseDirectionalCosts.vertCosts = vertCost;
pairwiseDirectionalCosts.horCosts = horCost;
if exist('mainDiagCost', 'var')
    pairwiseDirectionalCosts.mainDiagCosts = mainDiagCost;
end
if exist('secondDiagCost', 'var')
    pairwiseDirectionalCosts.secondDiagCosts = secondDiagCost;
end



end

