% example of usage of alphaExpansionRobustHighOrderPottsMex.m
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

dataTerm = [-10 10 0; 10 -10 0];
pairwiseTerm = sparse([0 1 2; 0 0 5; 0 0 0]);
ho = { [1; 2] };
hop = [ 10, 10 ];

[labels, energy,  energyPlot, timePlot] = alphaExpansionRobustHighOrderPottsMex(dataTerm, pairwiseTerm, ho, hop);


% % correct answer: 
% energy = -16; 
% labels = [1; 2; 2];

if ~isequal(energy, -16)
    warning('Wrong value of energy!')
end
if ~isequal(labels, [1; 2; 2])
    warning('Wrong value of labels!')
end

