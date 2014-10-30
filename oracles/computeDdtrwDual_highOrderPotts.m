function [dualValue, subgradient, primalLabeling] = computeDdtrwDual_highOrderPotts(dataCost, pairwiseCosts, dualVars, hoIds, hoP)
% computeDdtrwDual_highOrderPotts computes the value of the dual functon of the DD-TRW on a Potts Grid with robust P^n-Potts potentials
%
% The function minimizes the Lagrangian over binary variables Y given duals variables D:
% L(Y, D) = E_{vert}(Y^1) + E_{hor}(Y^2) +  E_{diag1}(Y^3) + E_{diag2}(Y^4) + \sum_{c} E_c(Y^c)
%   + \sum_i \sum_k d_{hor ik} (y^1_ik - y^2_ik) 
%   + \sum_i \sum_k d_{diag1 ik} (y^1_ik - y^3_ik)
%   + \sum_i \sum_k d_{diag2 ik} (y^1_ik - y^4_ik)
%   + \sum_{c} \sum_{i \in c} \sum_k d_{c, ik} ( y^1_{ik} - y^c_{ik} )
%
% [dualValue, subgradient, primalLabeling] = computeDdtrwDual_highOrderPotts(dataCost, pairwiseCosts, dualVars, hoIds, hoP)
%
% INPUT
%   dataCost   - unary potentials ( double[ numLabels x numNodes ])
%   pairwiseCosts - structure that contains Potts coefficients for the pairwise potentials. Fields:
%       vertCosts - vertical pairwise potentials ( double[ (heigth - 1), width ] ), heigth * width = numNodes
%       horCosts - horizontal pairwise potentials ( double[ heigth, (width - 1) ] ), 
%       mainDiagCosts - '\' diagonal potentials ( double[ (heigth - 1), (width - 1) ] ), 
%       secondDiagCosts - '/' diagonal potentials ( double[ (heigth - 1), (width - 1) ] ).
%   dualVars   - vector of dual varuables ( double[ (numNodes * numLabels *(numDirections - 1) + numLabels * totalSizeHo)  x 1 ])
% 	hoIds       - groups of edges, showing high-order potentials (cell[numHO, 1], each element - vector of indices)
% 	hoP         - parameters of Robust high-order potentials (double[numHO, 2]), each row gives \gamma_max and Q, \gamma =0; \gamma_max >= 0; Q >= 0;
%
% OUTPUT
%   dualValue - the value of the dual function
%   subgradient - value of subgradient
%   primalLabeling - the estimate of primal labeling
%
% Anton Osokin (firstname.lastname@gmail.com),  16.05.2013

if ~isnumeric(dataCost) || ~ismatrix(dataCost)
    error('computeDdtrwDual_highOrderPotts:badDataCost', 'dataCost should be a matrix  numLabels x numNodes');
end
dataCost = double(dataCost);
numNodes = size(dataCost, 2);
numLabels = size(dataCost, 1);

if ~isstruct(pairwiseCosts) || ~isfield(pairwiseCosts, 'vertCosts') || ~isfield(pairwiseCosts, 'horCosts')
    error('computeDdtrwDual_highOrderPotts:badPairwiseCost', 'pairwiseCosts should be a structure with fields vertCosts, horCosts, and possibly mainDiagCosts, secondDiagCosts ');
end
vertCost = pairwiseCosts.vertCosts;
horCost = pairwiseCosts.horCosts;

if ~isnumeric(vertCost) || ~ismatrix(vertCost) || (size(vertCost, 1) + 1) * size(vertCost, 2) ~= numNodes || any( vertCost(:) < 0)
    error('computeDdtrwDual_highOrderPotts:badVertCost', 'pairwiseCosts.vertCosts should be a matrix with nonnegative elements of size (heigth - 1) x width,  heigth * width = numNodes ');
end
gridSize = nan(2, 1);
gridSize(1) = size(vertCost, 1) + 1;
gridSize(2) = size(vertCost, 2);

if ~isnumeric(horCost) || ~ismatrix(horCost) || (size(horCost, 2) + 1) * size(horCost, 1) ~= numNodes || any( horCost(:) < 0) || size(horCost, 2) + 1 ~= gridSize(2) || size(horCost, 1) ~= gridSize(1)
    error('computeDdtrwDual_highOrderPotts:badHorCost', 'pairwiseCosts.horCosts should be a matrix with nonnegative elements of size heigth x (width - 1),  heigth * width = numNodes ');
end
numDualVars = numLabels * numNodes;

if isfield(pairwiseCosts, 'mainDiagCosts')
    mainDiagCost = pairwiseCosts.mainDiagCosts;
    if ~isnumeric(mainDiagCost) || ~ismatrix(mainDiagCost) || (size(mainDiagCost, 2) + 1 ~= gridSize(2)) || (size(mainDiagCost, 1) + 1 ~= gridSize(1)) || any( mainDiagCost(:) < 0)
        error('computeDdtrwDual_highOrderPotts:badDiagCost', 'pairwiseCosts.mainDiagCosts should be a matrix with nonnegative elements of size (heigth - 1) x (width - 1),  heigth * width = numNodes ');
    end
    
    numDualVars = numDualVars + numLabels * numNodes;
end

if isfield(pairwiseCosts, 'secondDiagCosts')
    secondDiagCost = pairwiseCosts.secondDiagCosts;
    if ~isnumeric(secondDiagCost) || ~ismatrix(secondDiagCost) || (size(secondDiagCost, 2) + 1 ~= gridSize(2)) || (size(secondDiagCost, 1) + 1 ~= gridSize(1)) || any( secondDiagCost(:) < 0)
        error('computeDdtrwDual_highOrderPotts:badDiagCost', 'pairwiseCosts.secondDiagCosts should be a matrix with nonnegative elements of size (heigth - 1) x (width - 1),  heigth * width = numNodes ');
    end
    
    numDualVars = numDualVars + numLabels * numNodes;
end

if ~iscell(hoIds) || ~iscolumn(hoIds) 
    error('computeDdtrwDual_highOrderPotts:badHoIds', 'hoIds should be a cell array of size numHighOrder x 1');
end
numHo = length(hoIds);

for iHo = 1 : numHo
    if ~isnumeric(hoIds{iHo}) || ~iscolumn(hoIds{iHo}) || any(hoIds{iHo}(:) < 1) || any(hoIds{iHo}(:) > numNodes) || any(round(hoIds{iHo}(:)) - hoIds{iHo}(:) > 1e-5)
        error('computeDdtrwDual_highOrderPotts:badHoIdsElement', ['Problem with ho potential ', num2str(iHo) ,'. Each cell of hoIds should be a vector of node indices incident to the current potential: integer numbers between 1 and numNodes']);
    end
    
    numDualVars = numDualVars + length( hoIds{iHo} ) * numLabels;
end

if ~isnumeric(hoP) || ~ismatrix(hoP) || size(hoP, 1) ~= numHo || size(hoP, 2) ~= 2 || any(hoP(:)) < 0
    error('computeDdtrwDual_highOrderPotts:badHoP', 'hoP should be a matrix numHO x 2, all elements should be positive, ');
end

if ~isnumeric(dualVars) || ~iscolumn(dualVars) || length(dualVars) ~= numDualVars
    error('computeDdtrwDual_highOrderPotts:badDualVars', ['dualVars should be a column vector of length compatible with the data: expected - ', num2str(numDualVars), ', obtained - ', num2str(length(dualVars)) ]);
end
dualVars = double(dualVars);


horVars = reshape(dualVars(1 : numLabels * numNodes), [numLabels, numNodes]);
% find other vars
shift = numLabels * numNodes;
hoVars = cell(numHo, 1);
for iHo = 1 : numHo
    hoVars{iHo} = dualVars(shift + 1 : shift  + numLabels * length( hoIds{iHo} ));
    shift = shift + numLabels * length( hoIds{iHo} );
    
    hoVars{iHo} = reshape( hoVars{iHo}, [numLabels, length( hoIds{iHo} )] );
end

vertVars = -horVars;
for iHo = 1 : numHo
    vertVars(:, hoIds{iHo}) = vertVars(:, hoIds{iHo}) - hoVars{iHo};
end

if exist('mainDiagCost', 'var')
    mainDiagVars = reshape( dualVars(shift + 1 : shift + numLabels * numNodes), [numLabels, numNodes]);
    shift = shift + numLabels * numNodes;
    vertVars = vertVars - mainDiagVars;
end

if exist('secondDiagCost', 'var')
    secondDiagVars = reshape( dualVars(shift + 1 : shift + numLabels * numNodes), [numLabels, numNodes]);
    shift = shift + numLabels * numNodes;
    vertVars = vertVars - secondDiagVars;
end


% dualVars = reshape(dualVars, [numLabels, numNodes]);

ids = reshape( 1 : prod(gridSize), [gridSize(1), gridSize(2)] );

tmp = 0.5 * dataCost;
dataCostHor = tmp + horVars;
dataCostVert = tmp + vertVars;

%% minimize horizontal chains
labelsHor = nan( gridSize(2), gridSize(1) );
energyHor = nan( gridSize(1), 1);
for iRow = 1 : gridSize(1)
    curUnary = dataCostHor(:, ids(iRow, :))';
    [energyHor(iRow), labelsHor(:, iRow)] = viterbiPottsMex(curUnary, horCost(iRow, :)');
end
labelsHor = labelsHor';

%% minimize vertical chains
labelsVert = nan( gridSize(1), gridSize(2) );
energyVert = nan( gridSize(2), 1);
for iCol = 1 : gridSize(2)
    curUnary = dataCostVert(:, ids(:, iCol))';
    [energyVert(iCol), labelsVert(:, iCol)] = viterbiPottsMex(curUnary, vertCost(:, iCol));
end

%% minimize main diag chains
if exist('mainDiagCost', 'var')
    labelsMainDiag = nan( gridSize(1), gridSize(2) );
    energyMainDiag = nan(  gridSize(1) + gridSize(2) - 1, 1);

    diagIds = toeplitz(1 : gridSize(1), [1, gridSize(1) + 1 : 1 : gridSize(1) + gridSize(2) - 1 ]);
    for iDiag = 1 : gridSize(1) + gridSize(2) - 1
        curIds = ids( diagIds == iDiag );
        curUnary = mainDiagVars(:, curIds)';
        
        curPairwiseChain = mainDiagCost( diagIds(1 : end - 1, 1 : end - 1) == iDiag );
        
        if ~isempty(curPairwiseChain)
            [energyMainDiag(iDiag), curLabels] = viterbiPottsMex(curUnary, curPairwiseChain);
            labelsMainDiag(diagIds == iDiag) = curLabels;
        else
            if length(curIds) > 1
                error('Something wrong has happened with the main diag potentials');
            end
            [energyMainDiag(iDiag), curLabels] = min(curUnary, [], 2);
            labelsMainDiag(curIds) = curLabels;
        end
    end
end

%% minimize second diag chains
if exist('secondDiagCost', 'var')
    labelsSecondDiag = nan( gridSize(1), gridSize(2) );
    energySecondDiag = nan(  gridSize(1) + gridSize(2) - 1, 1);

    diagIds = rot90(toeplitz(1 : gridSize(2), [1, gridSize(2) + 1 : 1 : gridSize(2) + gridSize(1) - 1 ]));
    for iDiag = 1 : gridSize(1) + gridSize(2) - 1
        curIds = ids( diagIds == iDiag );
        curUnary = secondDiagVars(:, curIds)';
        
        curPairwiseChain = secondDiagCost( diagIds(2 : end, 1 : end - 1) == iDiag );
        
        if ~isempty(curPairwiseChain)
            [energySecondDiag(iDiag), curLabels] = viterbiPottsMex(curUnary, curPairwiseChain);
            labelsSecondDiag(diagIds == iDiag) = curLabels;
        else
            if length(curIds) > 1
                error('Something wrong has happened with the second diag potentials');
            end
            [energySecondDiag(iDiag), curLabels] = min(curUnary, [], 2);
            labelsSecondDiag(curIds) = curLabels;
        end
    end
end


%% minimize high-order potentials
energyHo = nan( numHo, 1 );
labelsHo = cell( numHo, 1 );
for iHo = 1 : numHo
    curUnary = hoVars{iHo};
    
    [minUnaryVal, unaryLabels] = min( curUnary, [], 1 );
    
    curSize = length( hoIds{iHo} );
    
    gamma = hoP(iHo, 1);
    Q = hoP(iHo, 2);
    theta = gamma / Q;
    
    
    curBestEnergy = sum( minUnaryVal ) + gamma;
    curBestLabeling = unaryLabels(:);
    
    for iLabel = 1 : numLabels
        
        tmp = curUnary(iLabel, :) - minUnaryVal;
        
        [minDiff, id] = sort(tmp, 'ascend');
        
        [minVal, minPos] = min( [0, cumsum( minDiff )] + (curSize : -1 : 0) * theta );
        
        if (minVal < gamma) && (sum( minUnaryVal ) + minVal < curBestEnergy )
            % accept the labeling
            curBestEnergy = sum( minUnaryVal ) + minVal;
            curBestLabeling = iLabel * ones( curSize, 1 );
            curBestLabeling( id( minPos : 1 : end) ) = unaryLabels( id( minPos : 1 : end ) );
        end
    end
    
    labelsHo{ iHo } = curBestLabeling;
    energyHo( iHo ) = curBestEnergy;
end


%% compute results
dualValue = sum(energyVert) + sum(energyHor) + sum(energyHo);
if exist('mainDiagCost', 'var')
    dualValue = dualValue + sum(energyMainDiag);
end
if exist('secondDiagCost', 'var')
    dualValue = dualValue + sum(energySecondDiag);
end
primalLabeling = labelsVert(:);

% this appears to be very slow
% subgradient = zeros(numLabels, numNodes);
% for iLabel = 1 : numLabels
%     subgradient(iLabel, :) = subgradient(iLabel, :) + (labelsVert(:) == iLabel)' - (labelsHor(:)' == iLabel);
% end

subgradientHor = zeros(numLabels, numNodes);
curDiff = find( labelsVert ~= labelsHor );
for iNodeId = 1 : length(curDiff)
    iNode = curDiff(iNodeId);
    subgradientHor( labelsVert(iNode), iNode ) = -1;
    subgradientHor( labelsHor(iNode), iNode ) = 1;
end

subgradientHo = cell( numHo, 1 );
for iHo = 1 : numHo
    subgradientHo{ iHo } = zeros( numLabels, length( hoIds{iHo} ));
    
    curDiff = find( labelsHo{ iHo } ~= labelsVert( hoIds{iHo} ) );
    
    for iNodeId = 1 : length(curDiff)
        iNode = curDiff(iNodeId);
        subgradientHo{ iHo }( labelsVert( hoIds{iHo}(iNode) ), iNode ) = -1;
        subgradientHo{ iHo }( labelsHo{ iHo }(iNode), iNode ) = 1;
    end
    
    subgradientHo{ iHo } = subgradientHo{ iHo }(:);
end

subgradient = [subgradientHor(:); cat(1, subgradientHo{:})];

if exist('mainDiagCost', 'var')
    subgradientMainDiag = zeros(numLabels, numNodes);
    curDiff = find( labelsVert ~= labelsMainDiag );
    for iNodeId = 1 : length(curDiff)
        iNode = curDiff(iNodeId);
        subgradientMainDiag( labelsVert(iNode), iNode ) = -1;
        subgradientMainDiag( labelsMainDiag(iNode), iNode ) = 1;
    end
    subgradient = [subgradient; subgradientMainDiag(:)];
end


if exist('secondDiagCost', 'var')
    subgradientSecondDiag = zeros(numLabels, numNodes);
    curDiff = find( labelsVert ~= labelsSecondDiag );
    for iNodeId = 1 : length(curDiff)
        iNode = curDiff(iNodeId);
        subgradientSecondDiag( labelsVert(iNode), iNode ) = -1;
        subgradientSecondDiag( labelsSecondDiag(iNode), iNode ) = 1;
    end
    subgradient = [subgradient; subgradientSecondDiag(:)];
end

end
