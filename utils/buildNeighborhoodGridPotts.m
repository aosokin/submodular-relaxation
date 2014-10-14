function neighbors = buildNeighborhoodGridPotts(nRows, nCols, pottsWeight, connectivityType)
%buildNeighborhoodGridPotts creates a 4-connected or 8-connected graph with Potts weights
%
% neighbors = buildNeighborhoodGridPotts(nRows, nCols, pottsWeight, connectivityType)
%
% INPUT
%   nRows   - number of rows (double)
%   nCols   - number of columns (double)
%   pottsWeight - there are two options:
%                       1) potts weight of all edges (double)
%                       2) structure with fields vertCost, horCost and optionally diag1Cost and diag2Cost
%   connectivityType   -  string '4' or '8', default: '4'
%
% OUTPUT
%   neighbors  - sparse array of pairwise weights ( sparse double[ numNodes x numNodes ], numNodes = nRows * nCols )
%
% Anton Osokin (firstname.lastname@gmail.com),  25.04.2013

%% check data
if nargin == 3
    connectivityType = '4';
end

if ~ischar(connectivityType) || ~( strcmpi( connectivityType, '4') || strcmpi( connectivityType, '8') )
    error('buildNeighborhoodGridPotts:badConnectivityType', 'Connectivity type has to be string: "4" or "8"');
end

if ~isnumeric( nRows ) || ~isscalar(nRows) || nRows <= 0
    error('buildNeighborhoodGridPotts:badNRows', 'nRows should be an integer positive number');
end
nRows = double(nRows);
if abs( nRows - round(nRows)) > 1e-5
    error('buildNeighborhoodGridPotts:badNRowsNotInteger', 'nRows should be an integer positive number');
end

if ~isnumeric( nCols ) || ~isscalar(nCols) || nCols <= 0
    error('buildNeighborhoodGridPotts:badNCols', 'nCols should be an integer positive number');
end
nCols = double(nCols);
if abs( nCols - round(nCols)) > 1e-5
    error('buildNeighborhoodGridPotts:badNColsNotInteger', 'nCols should be an integer positive number');
end

if ~isnumeric( pottsWeight ) || ~isscalar(pottsWeight)
    if ~isstruct( pottsWeight ) || ~isfield(pottsWeight, 'vertCost') || ~isfield(pottsWeight, 'horCost')
        error('buildNeighborhoodGridPotts:badPottsWeight', 'pottsWeight should be a number or a specific structure');
    else
        if strcmpi( connectivityType, '8') && (~isfield(pottsWeight, 'diag1Cost') || ~isfield(pottsWeight, 'diag2Cost'))
            error('buildNeighborhoodGridPotts:wrongConnectivityTypeInPottsWeight', 'fields in pottsWeight should be consistent with connectivityType');
        end
        if  ~isnumeric( pottsWeight.vertCost ) || ~ismatrix(pottsWeight.vertCost) || size(pottsWeight.vertCost, 1) ~= nRows - 1 || size(pottsWeight.vertCost, 2) ~= nCols
            error('buildNeighborhoodGridPotts:badVertCost', 'field vertCost in pottsWeight should be of proper type and size');
        end
        if  ~isnumeric( pottsWeight.horCost ) || ~ismatrix(pottsWeight.horCost) || size(pottsWeight.horCost, 1) ~= nRows || size(pottsWeight.horCost, 2) ~= nCols - 1
            error('buildNeighborhoodGridPotts:badHorCost', 'field horCost in pottsWeight should be of proper type and size');
        end
        if strcmpi( connectivityType, '8')
            if  ~isnumeric( pottsWeight.diag1Cost ) || ~ismatrix(pottsWeight.diag1Cost) || size(pottsWeight.diag1Cost, 1) ~= nRows - 1 || size(pottsWeight.diag1Cost, 2) ~= nCols - 1
                error('buildNeighborhoodGridPotts:badDiag1Cost', 'field diag1Cost in pottsWeight should be of proper type and size');
            end
            if  ~isnumeric( pottsWeight.diag2Cost ) || ~ismatrix(pottsWeight.diag2Cost) || size(pottsWeight.diag2Cost, 1) ~= nRows - 1 || size(pottsWeight.diag2Cost, 2) ~= nCols - 1
                error('buildNeighborhoodGridPotts:badDiag2Cost', 'field diag2Cost in pottsWeight should be of proper type and size');
            end
        end
    end
else
    pottsWeightVal = double( pottsWeight );
    pottsWeight = struct;
    pottsWeight.vertCost = pottsWeightVal * ones(nRows - 1, nCols);
    pottsWeight.horCost = pottsWeightVal * ones(nRows, nCols - 1);
    if strcmpi( connectivityType, '8')
        pottsWeight.diag1Cost = pottsWeightVal * ones(nRows - 1, nCols - 1);
        pottsWeight.diag2Cost = pottsWeightVal * ones(nRows - 1, nCols - 1);
    end
end

if nRows == 1 && nCols == 1
    neighbors = sparse(0);
    return;
end

%% work
id = reshape( 1 : nRows * nCols, [nRows, nCols]);

% vertical
startEdge = id(1 : end - 1, :);
endEdge = id(2 : end, :);
startEdgeId = startEdge(:);
endEdgeId = endEdge(:);
edgeWeights = pottsWeight.vertCost(:);

% horizontal
startEdge = id(:, 1 : end - 1);
endEdge = id(:, 2 : end);
startEdgeId = [startEdgeId; startEdge(:)];
endEdgeId = [endEdgeId; endEdge(:)];
edgeWeights = [edgeWeights; pottsWeight.horCost(:)];

if strcmpi(connectivityType, '8');
    % diagonal \
    startEdge = id(1 : end - 1, 1 : end - 1);
    endEdge = id(2 : end, 2 : end);
    startEdgeId = [startEdgeId; startEdge(:)];
    endEdgeId = [endEdgeId; endEdge(:)];
    edgeWeights = [edgeWeights; pottsWeight.diag1Cost(:)];
    
    % diagonal /
    startEdge = id(1 : end - 1, 2 : end);
    endEdge = id(2 : end, 1 : end - 1);
    startEdgeId = [startEdgeId; startEdge(:)];
    endEdgeId = [endEdgeId; endEdge(:)];
    edgeWeights = [edgeWeights; pottsWeight.diag2Cost(:)];
end

% create symmentric sparse matrix
neighbors = sparse( [startEdgeId; endEdgeId], [endEdgeId; startEdgeId], [edgeWeights; edgeWeights] );

end
