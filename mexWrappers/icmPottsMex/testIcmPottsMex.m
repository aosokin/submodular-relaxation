function test_suite = testIcmPottsMex()
initTestSuite;
end

function testIcmPottsMex_badUnary
dataCost = struct;
neighbors = sparse([0, 1; 1, 0]);

f = @() icmPottsMex(dataCost, neighbors);
assertExceptionThrown(f, 'icmPottsMex:unaryPotentialsWrongType');
end


function testIcmPottsMex_badPairwise
dataCost = [1 0 0; 0 0 0];
neighbors = sparse([0, 1; 1, 0]);

f = @() icmPottsMex(dataCost, neighbors);
assertExceptionThrown(f, 'icmPottsMex:pairwisePotentialsWrongSize');
end

function testIcmPottsMex_badInitLabels
dataCost = [1 0 0; 0 0 0];
neighbors = sparse(eye(3));
labels = zeros( 3, 1 );

f = @() icmPottsMex(dataCost, neighbors, labels);
assertExceptionThrown(f, 'icmPottsMex:initLabelsWrongLabel');
end

function testIcmPottsMex_badMexIter
dataCost = [1 0 0; 0 0 0];
neighbors = sparse(eye(3));
labels = ones( 3, 1 );
maxIter = -1;

f = @() icmPottsMex(dataCost, neighbors, labels, maxIter);
assertExceptionThrown(f, 'icmPottsMex:maxNumIterNegativeValue');
end


function testComputeEnergyPotts_runIcmPottsMexZeroIter
dataCost = [ 0 1; 0 1; 0 1];
neighbors = sparse([0, 1; 0, 0]);
labels = [1; 2];
maxIter = 0;

[energy, newLabels] = icmPottsMex(dataCost, neighbors, labels, maxIter);

assertEqual(energy, 2);
assertEqual(newLabels, labels);
end

function testComputeEnergyPotts_runIcmPottsMexOneIter
dataCost = [ 0 1; 0 1; 0 1];
neighbors = sparse([0, 1; 0, 0]);
labels = [1; 2];
maxIter = 1;

[energy, newLabels] = icmPottsMex(dataCost, neighbors, labels, maxIter);

assertEqual(energy, 1);
assertEqual(newLabels, [2; 2]);
end

% function testComputeEnergyPotts_badDataCost
% dataCost = struct;
% neighbors = sparse([0, 1; 1, 0]);
% labels = [1; 2];
% 
% f = @() computeEnergyPotts(dataCost,neighbors, labels);
% assertExceptionThrown(f, 'computeEnergyPotts:badDataCost');
% end
% 
% function testComputeEnergyPotts_badNeighbors
% dataCost = [ 0 1; 0 1; 0 1];
% neighbors = sparse(eye(10));
% labels = [1; 2];
% 
% f = @() computeEnergyPotts(dataCost,neighbors, labels);
% assertExceptionThrown(f, 'computeEnergyPotts:badNeighbors');
% end
% 
% function testComputeEnergyPotts_badLabels
% dataCost = [ 0 1; 0 1; 0 1];
% neighbors = sparse([0, 1; 1, 0]);
% labels = [0; 1];
% 
% f = @() computeEnergyPotts(dataCost,neighbors, labels);
% assertExceptionThrown(f, 'computeEnergyPotts:badLabels');
% end
