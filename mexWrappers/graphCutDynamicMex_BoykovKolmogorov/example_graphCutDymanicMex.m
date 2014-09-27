dataTerms = [ 0 0.1; 0 0.1; 0 0.1];
pairwiseTerms = [1 2 0.5 0.7; ... % usual edge
                 1 3 2 -1; ... % one side reparametrization
                 3 2 -1 2]; % second side reparametrization

[energy, labels, graphHandle] = graphCutDynamicMex(dataTerms, pairwiseTerms);

if ~isequal(energy, -1.9)
    warning('Wrong value of energy!')
end
if ~isequal(labels, [1; 1; 0])
    warning('Wrong energy of labels!')
end

unaryUpdate = [3, 0, -1];

[energy, labels] = updateUnaryGraphCutDynamicMex(graphHandle, unaryUpdate);

if ~isequal(energy, -2.9)
    warning('Wrong value of energy!')
end
if ~isequal(labels, [1; 1; 0])
    warning('Wrong energy of labels!')
end

unaryUpdate = [1, 10, 0; 2 10 0; 3 -10 0];

[energy, labels] = updateUnaryGraphCutDynamicMex(graphHandle, unaryUpdate);

if ~isequal( round(energy * (1e+3)) / (1e+3), -5.8)
    warning('Wrong value of energy!')
end
if ~isequal(labels, [0; 0; 1])
    warning('Wrong energy of labels!')
end

deleteGraphCutDynamicMex( graphHandle );



