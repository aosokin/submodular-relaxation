function test_suite = testMrfMinimizeMex_time()
initTestSuite;
end


function testMrfMinimizeMex_time_badUnaryTerms
dataTerm = struct;
pairwiseTerm = [1 2 0 0];

f = @() mrfMinimizeMex_time(dataTerm, pairwiseTerm);
assertExceptionThrown(f, 'mrfMinimizeMex_time:error');
end

function testMrfMinimizeMex_time_badPairwiseTerms
dataTerm = eye(2);
pairwiseTerm = struct;

f = @() mrfMinimizeMex_time(dataTerm, pairwiseTerm);
assertExceptionThrown(f, 'mrfMinimizeMex_time:error');
end

function testMrfMinimizeMex_time_badMetric
dataTerm = eye(2);
pairwiseTerm = sparse(eye(2));
metric = eye(3);

f = @() mrfMinimizeMex_time(dataTerm, pairwiseTerm, metric);
assertExceptionThrown(f, 'mrfMinimizeMex_time:error');
end

function testMrfMinimizeMex_time_badOptions
dataTerm = eye(2);
pairwiseTerm = sparse(eye(2));
options = {1};

f = @() mrfMinimizeMex_time(dataTerm, pairwiseTerm, [], options);
assertExceptionThrown(f, 'mrfMinimizeMex_time:error');
end

function testMrfMinimizeMex_time_badOptionsMethod
dataTerm = eye(2);
pairwiseTerm = sparse(eye(2));
options = struct;
options.method = 1;

f = @() mrfMinimizeMex_time(dataTerm, pairwiseTerm, [], options);
assertExceptionThrown(f, 'mrfMinimizeMex_time:error');
end

function testMrfMinimizeMex_time_badOptionsMaxIter
dataTerm = eye(2);
pairwiseTerm = sparse(eye(2));
options = struct;
options.maxIter = -1;

f = @() mrfMinimizeMex_time(dataTerm, pairwiseTerm, [], options);
assertExceptionThrown(f, 'mrfMinimizeMex_time:error');
end

function testMrfMinimizeMex_time_badOptionsFuncEps
dataTerm = eye(2);
pairwiseTerm = sparse(eye(2));
options = struct;
options.funcEps = struct;

f = @() mrfMinimizeMex_time(dataTerm, pairwiseTerm, [], options);
assertExceptionThrown(f, 'mrfMinimizeMex_time:error');
end

function testMrfMinimizeMex_time_badOptionsVerbosity
dataTerm = eye(2);
pairwiseTerm = sparse(eye(2));
options = struct;
options.verbosity = -1;

f = @() mrfMinimizeMex_time(dataTerm, pairwiseTerm, [], options);
assertExceptionThrown(f, 'mrfMinimizeMex_time:error');
end

function testMrfMinimizeMex_time_badOptionsPrintMinIter
dataTerm = eye(2);
pairwiseTerm = sparse(eye(2));
options = struct;
options.printMinIter = -1;

f = @() mrfMinimizeMex_time(dataTerm, pairwiseTerm, [], options);
assertExceptionThrown(f, 'mrfMinimizeMex_time:error');
end

function testMrfMinimizeMex_time_badOptionsPrintIter
dataTerm = eye(2);
pairwiseTerm = sparse(eye(2));
options = struct;
options.printIter = -1;

f = @() mrfMinimizeMex_time(dataTerm, pairwiseTerm, [], options);
assertExceptionThrown(f, 'mrfMinimizeMex_time:error');
end

function testMrfMinimizeMex_time_runTrwsPottsMode
dataTerm = [-10 10 0; 10 -10 0];
pairwiseTerm = sparse([0 1 2; 0 0 5; 0 0 0]);
options = struct;
options.method = 'trw-s';

[S, E, LB, lbPlot, energyPlot, timePlot] = mrfMinimizeMex_time(dataTerm, pairwiseTerm, [], options);

assertEqual(S, [1; 2; 2]);
assertEqual(E, -17); 
assertEqual(LB, -17); 
assertTrue(min(energyPlot) == -17);
assertTrue(max(lbPlot) == -17);
assertTrue( all(diff(timePlot) <= 0) );
end

function testMrfMinimizeMex_time_runTrwsGeneralMode
dataTerm = [-10 10 0; 10 -10 0];
pairwiseTerm = sparse([0 1 2; 0 0 5; 0 0 0]);
options = struct;
options.method = 'trw-s';
metric = ones(2) - eye(2);

[S, E, LB, lbPlot, energyPlot, timePlot] = mrfMinimizeMex_time(dataTerm, pairwiseTerm, metric, options);

assertEqual(S, [1; 2; 2]);
assertEqual(E, -17); 
assertEqual(LB, -17); 
assertTrue(min(energyPlot) == -17);
assertTrue(max(lbPlot) == -17);
assertTrue( all(diff(timePlot) <= 0) );
end

function testMrfMinimizeMex_time_runBpPottsMode
dataTerm = [-10 10 0; 10 -10 0];
pairwiseTerm = sparse([0 1 2; 0 0 5; 0 0 0]);
options = struct;
options.method = 'bp';

[S, E, LB, lbPlot, energyPlot, timePlot] = mrfMinimizeMex_time(dataTerm, pairwiseTerm, [], options);

assertEqual(S, [1; 2; 2]);
assertEqual(E, -17); 
assertEqual(LB, NaN); 
assertTrue(min(energyPlot) == -17);
assertTrue(all(isnan(lbPlot)));
assertTrue( all(diff(timePlot) <= 0) );
end

function testMrfMinimizeMex_time_runBpGeneralMode
dataTerm = [-10 10 0; 10 -10 0];
pairwiseTerm = sparse([0 1 2; 0 0 5; 0 0 0]);
options = struct;
options.method = 'bp';
metric = ones(2) - eye(2);

[S, E, LB, lbPlot, energyPlot, timePlot] = mrfMinimizeMex_time(dataTerm, pairwiseTerm, metric, options);

assertEqual(S, [1; 2; 2]);
assertEqual(E, -17); 
assertEqual(LB, NaN); 
assertTrue(min(energyPlot) == -17);
assertTrue(all(isnan(lbPlot)));
assertTrue( all(diff(timePlot) <= 0) );
end

