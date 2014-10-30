function [energy, dataTerm, smoothnessTerm, highOrderTerm ] = computeEnergy_highOrderPotts(dataCost, neighbors, labels, hoIds, hoP)
%computeEnergy_highOrderPotts computes the value of energy with pairwise Potts potentials and robust high-order Potts potentials
%
% E(X) = \sum_i \sum_p U_{ip} [x_i == p]
%       + \sum_{ij} P_{ij} [x_i <> x_j]
%       + \sum_{c \in C} ( gamma_max + \sum_p  min{ (|c| - \sum_{i \in c } [x_i == p]) * gamma_max / Q - gamma_max, 0} )
%       + \sum_i d_i ( \sum_p y_{ip} - 1)
%
% [energy, dataTerm, smoothnessTerm] = computeEnergy_highOrderPotts(dataCost, neighbors, labels)
%
% INPUT
%   dataCost   - unary potentials ( double[ numLabels x numNodes ])
%   neighbors  - paiwise Potts potentials ( sparse double[ numNodes x numNodes ]). The function uses only upper triangle of this matrix.
%   labels - the labeling ( double[ numNodes x 1 ], 1-based indexing for labels)
% 	hoIds       - groups of nodes, showing high-order potentials (cell[numHO, 1], each element - vector of indices)
% 	hoP         - parameters of Robust high-order potentials (double[numHO, 2]), each row gives \gamma_max and Q; \gamma_max >= 0; Q >= 0;
%
% OUTPUT
%   energy - the energy value
%   dataTerm - the unary term
%   smoothnessTerm - the pairwise term
%   hoTerm - the high-order term
%
% Anton Osokin (firstname.lastname@gmail.com),  14.05.2013

if ~isnumeric(dataCost) || ~ismatrix(dataCost)
    error('computeEnergy_highOrderPotts:badDataCost', 'dataCost should be a matrix numLabels x numNodes');
end
dataCost = double(dataCost);
numLabels = size(dataCost, 1);
numNodes = size(dataCost, 2);

if ~isnumeric(neighbors) || ~ismatrix(neighbors) || ~issparse(neighbors) || size(neighbors, 1) ~= numNodes || size(neighbors, 2) ~= numNodes
    error('computeEnergy_highOrderPotts:badNeighbors', 'neighbors should be a sparse matrix numNodes x numNodes');
end

if ~iscell(hoIds) || ~iscolumn(hoIds) 
    error('computeEnergy_highOrderPotts:badHoIds', 'hoIds should be a cell array of size numHighOrder x 1');
end
numHO = length(hoIds);

for iHO = 1 : numHO
    if ~isnumeric(hoIds{iHO}) || ~iscolumn(hoIds{iHO}) || any(hoIds{iHO}(:) < 1) || any(hoIds{iHO}(:) > numNodes) || any(round(hoIds{iHO}(:)) - hoIds{iHO}(:) > 1e-5)
        error('computeEnergy_highOrderPotts:badHoIdsElement', ['Problem with ho potential ', num2str(iHO) ,'. Each cell of hoIds should be a vector of node indices incident to the current potential: integer numbers between 1 and numNodes']);
    end
end

if ~isnumeric(hoP) || ~ismatrix(hoP) || size(hoP, 1) ~= numHO || size(hoP, 2) ~= 2 || any(hoP(:)) < 0
    error('computeEnergy_highOrderPotts:badHoP', 'hoP should be a matrix numHO x 2, all elements should be positive, ');
end

if ~isnumeric(labels) || ~iscolumn(labels) || length(labels) ~= numNodes || any(~ismember( unique(labels), 1 : numLabels))
    error('computeEnergy_highOrderPotts:badLabels', 'labels should be a column vector of length numNodes, 1-indexing for label IDs');
end

dataTerm = sum(dataCost(labels(:) + (0 : numNodes - 1)' * numLabels));

[row,col,v] = find(neighbors);
newRow = row(row < col);
newCol = col(row < col);
newV = v(row < col);

smoothnessTerm = sum((labels(newRow) ~= labels(newCol)) .* newV);

highOrderTerm = 0;

numHO = length(hoIds);
for iHO = 1 : numHO
    curLabels = labels( hoIds{iHO} );
    curSize = length( curLabels );
    curStats = accumarray( curLabels, ones(size(curLabels)), [numLabels, 1]);  
    
    Q = hoP( iHO, 2);
    gamma = hoP( iHO, 1);
    theta = gamma / Q;
    
    highOrderTerm = highOrderTerm + min( min( (curSize - curStats) * theta ), gamma);
end

energy = dataTerm + smoothnessTerm + highOrderTerm ;

end

