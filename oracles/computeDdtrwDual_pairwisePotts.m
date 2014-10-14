function [dualValue, subgradient, primalLabeling] = computeDdtrwDual_pairwisePotts(dataCost, vertCost, horCost, dualVars)
%computeDdtrwDualPotts computes the value of the dual functon of the grid DD-TRW
%
% Please cite the following paper if you use this method:
%   N. Komodakis, N. Paragios, and G. Tziritas,
%   MRF energy minimization and beyond via dual decomposition.
%   IEEE TPAMI, vol. 33, no. 3, pp. 531–552, 2011.
%
% The function minimizes the Lagrangian over binary variables Y given duals variables D:
% L(Y, D) = E_{vert}(Y^1) + E_{hor}(Y^2) + \sum_i \sum_k D_{ik} ( Y^1_{ik} - Y^2_{ik} )
%
% [dualValue, subgradient, primalLabeling] = computeDdtrwDualPotts(dataCost, vertCost, horCost, dualVars)
%
% INPUT
%   dataCost   - unary potentials ( double[ numLabels x numNodes ])
%   vertCost - vertical pairwise potentials ( double[ (heigth - 1) * width ] ), heigth * width = numNodes 
%   horCost - horizontal pairwise potentials ( double[ heigth * (width - 1) ] ), heigth * width = numNodes
%   dualVars   - vector of dual varuables ( double[ numNodes*numLabels x 1 ])
%
% OUTPUT
%   dualValue - the value of the dual function
%   subgradient - value of subgradient
%   primalLabeling - the estimate of primal labeling
%
% Depends on mexWrappers/viterbiPottsMex
%
% To construct vertCost and horCost from sparse matrix neighbors use separateVertHorCosts.m
%
% Anton Osokin (firstname.lastname@gmail.com),  16.05.2013

if ~isnumeric(dataCost) || ~ismatrix(dataCost)
    error('computeDdtrwDual_pairwisePotts:badDataCost', 'dataCost should be a matrix  numLabels x numNodes');
end
dataCost = double(dataCost);
numNodes = size(dataCost, 2);
numLabels = size(dataCost, 1);

if ~isnumeric(vertCost) || ~ismatrix(vertCost) || (size(vertCost, 1) + 1) * size(vertCost, 2) ~= numNodes 
    error('computeDdtrwDual_pairwisePotts:badVertCost', 'vertCost a matrix with real-valued elements of size (heigth - 1) x width,  heigth * width = numNodes ');
end
gridSize = nan(2, 1);
gridSize(1) = size(vertCost, 1) + 1;
gridSize(2) = size(vertCost, 2);

if ~isnumeric(horCost) || ~ismatrix(horCost) || (size(horCost, 2) + 1) * size(horCost, 1) ~= numNodes 
    error('computeDdtrwDual_pairwisePotts:badHorCost', 'horCost a matrix with real-valued elements of size heigth x (width - 1),  heigth * width = numNodes ');
end
if size(horCost, 2) + 1 ~= gridSize(2) || size(horCost, 1) ~= gridSize(1)
    error('computeDdtrwDual_pairwisePotts:vertHorCostMismatch', 'vertCost and horCost are not compatible');
end


if ~isnumeric(dualVars) || ~iscolumn(dualVars) || length(dualVars) ~= numNodes * numLabels
    error('computeDdtrwDual_pairwisePotts:badDualVars', 'dualVars should be a column vector of length numNodes*numLabels');
end
dualVars = double(dualVars);
dualVars = reshape(dualVars, [numLabels, numNodes]);

ids = reshape( 1 : prod(gridSize), [gridSize(1), gridSize(2)] );

tmp = 0.5 *dataCost;
dataCostHor = tmp - dualVars;
dataCostVert = tmp + dualVars;

%% compute horizontal chains
labelsHor = nan( gridSize(2), gridSize(1) );
energyHor = nan( gridSize(1), 1);
for iRow = 1 : gridSize(1)
    curUnary = dataCostHor(:, ids(iRow, :))';
    [energyHor(iRow), labelsHor(:, iRow)] = viterbiPottsMex(curUnary, horCost(iRow, :)');
end
labelsHor = labelsHor';

%% compute vertical chains
labelsVert = nan( gridSize(1), gridSize(2) );
energyVert = nan( gridSize(2), 1);
for iCol = 1 : gridSize(2)
    curUnary = dataCostVert(:, ids(:, iCol))';
    [energyVert(iCol), labelsVert(:, iCol)] = viterbiPottsMex(curUnary, vertCost(:, iCol));
end

%% compute results
dualValue = sum(energyVert) + sum(energyHor);
primalLabeling = labelsVert(:);

% this appears to be very slow
% subgradient = zeros(numLabels, numNodes);
% for iLabel = 1 : numLabels
%     subgradient(iLabel, :) = subgradient(iLabel, :) + (labelsVert(:) == iLabel)' - (labelsHor(:)' == iLabel);
% end

subgradient = zeros(numLabels, numNodes);
curDiff = find( labelsVert ~= labelsHor );
for iNodeId = 1 : length(curDiff)
    iNode = curDiff(iNodeId);
    subgradient( labelsVert(iNode), iNode ) = 1;
    subgradient( labelsHor(iNode), iNode ) = -1;
end

subgradient = subgradient(:);

end
