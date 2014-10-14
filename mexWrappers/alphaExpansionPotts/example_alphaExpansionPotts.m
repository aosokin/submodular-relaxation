% example of usage of alphaExpansionPotts.m
%
% Anton Osokin (firstname.lastname@gmail.com),  25.04.2013

dataCost = [0 0 0 0 0; 1 -1 0 0 0];

neighbors = sparse([1; 1; 3; 3], [5; 3; 4; 5], [1; 10; 1; 1], 5, 5);

roundingPrecision = 1;

[labels, energy, time, energyPlot, timePlot] = alphaExpansionPotts(dataCost, neighbors, roundingPrecision);

% % correct answer: 
% energy = -1; 
% labels = [1; 2; 1; 1; 1];

if ~isequal(energy, -1)
    warning('Wrong value of energy!')
end
if ~isequal(labels, [1; 2; 1; 1; 1])
    warning('Wrong value of labels!')
end
