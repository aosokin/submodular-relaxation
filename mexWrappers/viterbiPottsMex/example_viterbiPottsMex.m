% example of usage of package viterbiPottsMex
%
% Anton Osokin (firstname.lastname@gmail.com),  22.05.2013

unary = [0 1 1; 1 0 1; 0 0 1; 0 0 -10];
costs = [2; 2; 3];

[energy, labels] = viterbiPottsMex(unary, costs);

% correct answer: energy = -7; labels = [3; 3; 3; 3];
if ~isequal(energy, -7)
    warning('Wrong value of energy!')
end
if ~isequal(labels, [3; 3; 3; 3])
    warning('Wrong value of labels!')
end
