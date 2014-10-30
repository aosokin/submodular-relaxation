function [dualValue, subgradient, primalLabeling] = computeSmrDualDynamic_highOrderPotts(dataCost, neighbors, dualVars, hoIds, hoP)
%computeSmrDualDynamic_highOrderPotts computes the value of the SMR dual function for energy with associative pairwise Potts potentials and robust high-order Potts potentials
%
% The function minimizes the Lagrangian w.r.t. binary variables  Y given duals variables D:
% L(Y, D) = \sum_i \sum_p U_{ip} y_{ip} 
%       + \sum_{ij} P_{ij} \sum_{p} 0.5 * ( [ y_{ip} == 1][y_{ip} == 0] + [ y_{ip} == 0][y_{ip} == 1] ) 
%       + \sum_{c \in C} ( gamma_max + \sum_p  min{ (|c| - \sum_{i \in c } y_{ip}) * gamma_max / Q - gamma_max, 0} )
%       + \sum_i d_i ( \sum_p y_{ip} - 1)
%
%   this function makes use of dynamic graph cuts to compute the updates faster
%   the following global variables are used: computeSmrDualDynamic_highOrderPotts_graphHandle, computeSmrDualDynamic_highOrderPotts_lastPoint
%           
%
% [dualValue, subgradient, primalLabeling]= computeSmrDualDynamic_highOrderPotts(dataCost, neighbors, dualVars)
%
% INPUT
%   dataCost   - unary potentials ( double[ numLabels x numNodes ])
%   neighbors  - paiwise Potts potentials ( sparse double[ numNodes x numNodes ]). The function uses only upper triangle of this matrix.
%   dualVars   - vector of dual varuables ( double[ numNodes x 1 ])
% 	hoIds       - groups of edges, showing high-order potentials (cell[numHO, 1], each element - vector of indices)
% 	hoP         - parameters of Robust high-order potentials (double[numHO, 2]), each row gives \gamma_max and Q; \gamma_max >= 0; Q >= 0;
%
% OUTPUT
%   dualValue - the value of the dual function
%   subgradient - value of subgradient
%   primalLabeling - the estimate of primal labeling
%
% Anton Osokin (firstname.lastname@gmail.com),  23.05.2013

% check for the mex-file
if ~isnumeric(dataCost) || ~ismatrix(dataCost)
    error('computeSmrDualDynamic_highOrderPotts:badDataCost', 'dataCost should be a matrix  numLabels x numNodes');
end
dataCost = double(dataCost);
numNodes = size(dataCost, 2);
numLabels = size(dataCost, 1);

if ~isnumeric(neighbors) || ~ismatrix(neighbors) || ~issparse(neighbors) || size(neighbors, 1) ~= numNodes || size(neighbors, 2) ~= numNodes
    error('computeSmrDualDynamic_highOrderPotts:badNeighbors', 'neighbors should be a sparse matrix numNodes x numNodes');
end

if ~isnumeric(dualVars) || ~iscolumn(dualVars) || length(dualVars) ~= numNodes
    error('computeSmrDualDynamic_highOrderPotts:badDualVars', 'dualVars should be a column vector of length numNodes');
end
dualVars = double(dualVars);

if ~iscell(hoIds) || ~iscolumn(hoIds) 
    error('computeSmrDualDynamic_highOrderPotts:badHoIds', 'hoIds should be a cell array of size numHighOrder x 1');
end
numHo = length(hoIds);

for iHO = 1 : numHo
    if ~isnumeric(hoIds{iHO}) || ~iscolumn(hoIds{iHO}) || any(hoIds{iHO}(:) < 1) || any(hoIds{iHO}(:) > numNodes) || any(round(hoIds{iHO}(:)) - hoIds{iHO}(:) > 1e-5)
        error('computeSmrDualDynamic_highOrderPotts:badHoIdsElement', ['Problem with ho potential ', num2str(iHO) ,'. Each cell of hoIds should be a vector of node indices incident to the current potential: integer numbers between 1 and numNodes']);
    end
end

if ~isnumeric(hoP) || ~ismatrix(hoP) || size(hoP, 1) ~= numHo || size(hoP, 2) ~= 2 || any(hoP(:)) < 0
    error('computeSmrDualDynamic_highOrderPotts:badHoP', 'hoP should be a matrix numHO x 2, all elements should be positive, ');
end

subEnergy = nan(numLabels, 1);
labelsQp = nan(numNodes, numLabels);

dynamicCutRebuildNumber = 50;

global computeSmrDualDynamic_highOrderPotts_graphHandle
global computeSmrDualDynamic_highOrderPotts_lastPoint
global computeSmrDualDynamic_highOrderPotts_dynamicNumber 

sumGamma = 0;
for iHo = 1 : numHo
    gamma_max = hoP(iHo, 1);
    sumGamma = sumGamma + gamma_max;
end

if isempty(computeSmrDualDynamic_highOrderPotts_graphHandle) || isempty(computeSmrDualDynamic_highOrderPotts_lastPoint) || isempty(computeSmrDualDynamic_highOrderPotts_dynamicNumber)...
        || ~iscell(computeSmrDualDynamic_highOrderPotts_graphHandle) || numel( computeSmrDualDynamic_highOrderPotts_graphHandle ) ~= numLabels ...
        || ~iscolumn(computeSmrDualDynamic_highOrderPotts_lastPoint) || length( computeSmrDualDynamic_highOrderPotts_lastPoint ) ~=  numNodes ...
        || ~isscalar(computeSmrDualDynamic_highOrderPotts_dynamicNumber) || ~isnumeric(computeSmrDualDynamic_highOrderPotts_dynamicNumber) ...
        || mod( computeSmrDualDynamic_highOrderPotts_dynamicNumber, dynamicCutRebuildNumber) == 0
    % remove all graphs if left
    if ~isempty(computeSmrDualDynamic_highOrderPotts_graphHandle) 
        for iLabel = 1 : length(computeSmrDualDynamic_highOrderPotts_graphHandle);
            deleteGraphCutDynamicMex( computeSmrDualDynamic_highOrderPotts_graphHandle{iLabel} )
        end
    end
    
    % construct edges for a graph cut
    [rowNeighbor, colNeighbor, weightNeighbor] = find(neighbors);
    deleteMask = rowNeighbor >= colNeighbor;
    rowNeighbor( deleteMask ) = [];
    colNeighbor( deleteMask ) = [];
    weightNeighbor( deleteMask ) = [];
    nonTermEdgesWeights = [rowNeighbor, colNeighbor, 0.5 * weightNeighbor, 0.5 * weightNeighbor];
    
    % add high-order cliques
    totalSize = 0;
    for iHo = 1 : numHo
        totalSize = totalSize + length( hoIds{iHo} );
    end
    extraUnary = nan(numHo, 2);
    extraPairwise = nan( totalSize, 4);
    
    shift = 0;
    sumGamma = 0;
    for iHo = 1 : numHo
        gamma_max = hoP(iHo, 1);
        Q = hoP(iHo, 2);
        
        sumGamma = sumGamma + gamma_max;

        extraUnary(iHo, 1) = -gamma_max;
        extraUnary(iHo, 2) = 0;
        
        extraPairwise( shift + 1 : shift + length(hoIds{iHo}), 1) = numNodes + iHo;
        extraPairwise( shift + 1 : shift + length(hoIds{iHo}), 2) = hoIds{iHo};
        extraPairwise( shift + 1 : shift + length(hoIds{iHo}), 3) = 0;
        extraPairwise( shift + 1 : shift + length(hoIds{iHo}), 4) = gamma_max / Q;
        
        shift = shift +  length(hoIds{iHo});
    end

    % construct unary terms for a graph cut
    termEdgeWeight  = dataCost';

    nonTermEdgesWeights = [nonTermEdgesWeights; extraPairwise];
    
    % store a point
    computeSmrDualDynamic_highOrderPotts_lastPoint = dualVars;
    computeSmrDualDynamic_highOrderPotts_graphHandle = cell(numLabels, 1);
    computeSmrDualDynamic_highOrderPotts_dynamicNumber = 1;

    for iLabel = 1 : numLabels
        
        curUnary = [termEdgeWeight(:, iLabel) + dualVars, zeros(numNodes, 1); extraUnary];
        
        [subEnergy(iLabel), curLabels ,...
            computeSmrDualDynamic_highOrderPotts_graphHandle{ iLabel }] = ...
            graphCutDynamicMex(curUnary, nonTermEdgesWeights);
        
        labelsQp(:, iLabel) = curLabels( 1 : numNodes);
    end
else
    pointDifference = dualVars - computeSmrDualDynamic_highOrderPotts_lastPoint;
    
    changeMask = pointDifference ~= 0;
    numChanges = sum( changeMask );
    
    unaryUpdate = [find(changeMask), pointDifference( changeMask ), zeros( numChanges, 1 )];
    for iLabel = 1 : numLabels
        
        [subEnergy(iLabel), curLabels] = updateUnaryGraphCutDynamicMex( computeSmrDualDynamic_highOrderPotts_graphHandle{ iLabel }, unaryUpdate );
        
        labelsQp(:, iLabel) = curLabels( 1 : numNodes);
    end
    computeSmrDualDynamic_highOrderPotts_lastPoint = dualVars;
    computeSmrDualDynamic_highOrderPotts_dynamicNumber = computeSmrDualDynamic_highOrderPotts_dynamicNumber + 1;
end

% compute the lower bound
dualValue = sum(subEnergy) - sum(dualVars) + sumGamma;

% get the primal estimate
if nargout > 2
    [~, primalLabeling] = max(labelsQp, [], 2);
end

%Compute subgradient
subgradient = sum(labelsQp, 2) - 1;

end
