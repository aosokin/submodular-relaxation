function [labels, energy, time, energyPlot, timePlot] = alphaExpansionPotts(dataCost, neighbors, roundingPrecision)
% alphaExpansionPotts runs alpha-expansion on an energy with Potts pairwise potentials
%   CAUTION: alphaExpansionPotts rounds both unary and pairwise potentials up to precision 1/roundingMultiplyer
%
% [labels, energy, time] = alphaExpansionPotts(dataCost, neighbors)
% [labels, energy, time, energyPlot, timePlot] = alphaExpansionPotts(dataCost, neighbors)
%
% INPUT
%   dataCost   - unary potentials ( double[ numLabels x numNodes ])
%   neighbors  - paiwise Potts potentials ( sparse double[ numNodes x numNodes ])
%   roundingPrecision - potentials are devided by roundingPrecision and then truncated ( double ), default: 1e-2
%
% OUTPUT
%   labels - labeling ( double[ numNodes, 1 ] )
%   energy - energy of labels ( double )
%   time - time ellapsed ( double )
%   timePlot, energyPlot - data for ploting energy against iterations ( double[ numIterations, 1 ] )
%                               plot( timePlot, energyPlot )    
%
% Depends on GCO library v 3.0: http://vision.csd.uwo.ca/code/gco-v3.0.zip
%
% Anton Osokin (firstname.lastname@gmail.com),  25.04.2013

if ~exist('dataCost', 'var')
    error('alphaExpansionPotts:noDataCost', 'dataCost input parameter is not specified');
end
if ~exist('neighbors', 'var')
    error('alphaExpansionPotts:noNeighbors', 'neighbors input parameter is not specified');
end
if nargin < 3
    roundingPrecision = 0.01;
end

% check the input
if ~isnumeric(dataCost) || ~ismatrix(dataCost)
    error('alphaExpansionPotts:badDataCost', 'dataCost should be a matrix  numLabels x numNodes');
end
dataCost = double(dataCost);
numNodes = size(dataCost, 2);
numLabels = size(dataCost, 1);

if ~isnumeric(neighbors) || ~ismatrix(neighbors) || ~issparse(neighbors) || size(neighbors, 1) ~= numNodes || size(neighbors, 2) ~= numNodes
    error('alphaExpansionPotts:badNeighbors', 'neighbors should be a sparse matrix numNodes x numNodes');
end
if ~isnumeric(roundingPrecision) || ~isscalar(roundingPrecision) || roundingPrecision <= 0
    error('alphaExpansionPotts:badRoundingPrecision', 'roundingPrecision should be a positive number');
end
roundingPrecision = double(roundingPrecision);

tStart = tic;
energyPlot = zeros(0, 1);
timePlot = zeros(0, 1);

Handle = GCO_Create(numNodes, numLabels);
GCO_SetDataCost(Handle, int32(dataCost / roundingPrecision));
GCO_SetNeighbors(Handle, round(neighbors / roundingPrecision));
GCO_SetSmoothCost(Handle, int32(ones(numLabels, numLabels) - eye(numLabels)));

if nargout > 3 
    % run alpha-expansion iteration by iteration
    energyPlot = nan(1000, 1);
    timePlot = nan(1000, 1);
    energyPlot(1) = double( GCO_ComputeEnergy(Handle) );
    timePlot(1) = toc(tStart);
    iter = 1;
    
    while iter == 1 || energyPlot(iter) < energyPlot(iter - 1)
        iter = iter + 1;
        energyPlot(iter) = double( GCO_Expansion(Handle, 1) );
        timePlot(iter) = toc(tStart);
    end 
    
    energyPlot( isnan(energyPlot) ) = [];
    timePlot( isnan(timePlot) ) = [];
    energy = min(energyPlot);
else
    energy = GCO_Expansion(Handle);
end

labels = double( GCO_GetLabeling(Handle) );
energy = double(energy);
GCO_Delete(Handle);

energyPlot = energyPlot * roundingPrecision;
energy = energy * roundingPrecision;

time = toc(tStart);

end
