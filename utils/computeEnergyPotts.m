function [energy, dataTerm, smoothnessTerm] = computeEnergyPotts(dataCost, neighbors, labels)
%computeEnergyPotts computes the value of pairwise energy with Potts potentials
% 
% [energy, dataTerm, smoothnessTerm] = computeEnergyPotts(dataCost, neighbors, labels)
%
% INPUT
%   dataCost   - unary potentials ( double[ numLabels x numNodes ])
%   neighbors  - paiwise Potts potentials ( sparse double[ numNodes x numNodes ]). The function uses only upper triangle of this matrix.
%   labels - the labeling ( double[ numNodes x 1 ], 1-based indexing for labels)
%
% OUTPUT
%   energy - the energy value
%   dataTerm - the unary term
%   smoothnessTerm - the pairwise term
%
% Anton Osokin (firstname.lastname@gmail.com),  14.05.2013

if ~isnumeric(dataCost) || ~ismatrix(dataCost)
    error('computeEnergyPotts:badDataCost', 'dataCost should be a matrix  numLabels x numNodes');
end
dataCost = double(dataCost);
numLabels = size(dataCost, 1);
numNodes = size(dataCost, 2);

if ~isnumeric(neighbors) || ~ismatrix(neighbors) || ~issparse(neighbors) || size(neighbors, 1) ~= numNodes || size(neighbors, 2) ~= numNodes
    error('computeEnergyPotts:badNeighbors', 'neighbors should be a sparse matrix numNodes x numNodes');
end

if ~isnumeric(labels) || ~iscolumn(labels) || length(labels) ~= numNodes || any(~ismember( unique(labels), 1 : numLabels))
    error('computeEnergyPotts:badLabels', 'labels should be a column vector of length numNodes, 1-indexing for label IDs');
end

dataTerm = sum(dataCost(labels(:) + (0 : numNodes - 1)' * numLabels));

[row,col,v] = find(neighbors);
newRow = row(row < col);
newCol = col(row < col);
newV = v(row < col);

smoothnessTerm = sum((labels(newRow) ~= labels(newCol)) .* newV);

energy = dataTerm + smoothnessTerm;

end

