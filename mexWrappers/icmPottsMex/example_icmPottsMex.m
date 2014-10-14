% example of usage of package icmPottsMex
%
% Anton Osokin (firstname.lastname@gmail.com),  14.10.2014

dataCost = [ 0 1; 0 1; 0 1];
neighbors = sparse([0, 1; 0, 0]);
labels = [1; 2];
maxIter = 1;

[energy, labels] = icmPottsMex(dataCost, neighbors, labels, maxIter);

% % correct answer: 
% energy = 1; 
% labels = [2; 2];

if ~isequal(energy, 1)
    warning('Wrong value of energy!')
end
if ~isequal(labels, [2; 2])
    warning('Wrong value of labels!')
end

