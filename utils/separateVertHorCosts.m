function [ vertCost, horCost, newNeighbors ] = separateVertHorCosts( neighbors, gridSize )
% separateVertHorCosts extracts cost for vertical and horizontal edges from sparse matrix
% 
% 	Anton Osokin (firstname.lastname@gmail.com),  14.10.2014

if ~isnumeric(gridSize) || ~isvector(gridSize) || length(gridSize) ~= 2 
    error('separateVertHorCosts:badGridSize', 'gridSize should be a column vector of length 2');
end
gridSize = double(gridSize);
numNodes = prod( gridSize );

if ~isnumeric(neighbors) || ~ismatrix(neighbors) || ~issparse(neighbors) || size(neighbors, 1) ~= numNodes || size(neighbors, 2) ~= numNodes
    error('separateVertHorCosts:badNeighbors', 'neighbors should be a sparse matrix numNodes x numNodes');
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

if sum( vertMask ) + sum( horMask ) ~= length(neighRow)
    warning('separateVertHorCosts:wrongNeighborhoodSystem', 'The graph is not 4-connected grid!!!');
end

vertCost = reshape(neighWeight(vertMask), [gridSize(1) - 1, gridSize(2)]);
horCost = reshape(neighWeight(horMask), [gridSize(1), gridSize(2) - 1]);

toKeepMask = vertMask | horMask;
newNeighbors = sparse(neighRow(toKeepMask), neighCol(toKeepMask), neighWeight(toKeepMask), numNodes, numNodes);


end

