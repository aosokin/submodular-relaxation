function tests = testViterbiPottsMex()
tests = functiontests(localfunctions);
end

function testViterbiPottsMex_badUnary(testCase)
unary = {};
costs = (1 : 10)';

f = @() viterbiPottsMex(unary, costs);
verifyError(testCase, f, 'viterbiPottsMex:unaryPotentialsWrongType');
end

function testViterbiPottsMex_badCosts(testCase)
unary = ones(9, 5);
unary(:, 1) = -1;
costs = -10 * ones(8, 1);

[energy, labels] = viterbiPottsMex(unary, costs);
verifyEqual(testCase, energy, -81);
verifyEqual(testCase, labels, [1; 2; 1; 2; 1; 2; 1; 2; 1]);
end

function testViterbiPottsMex_normalRun(testCase)
unary = [0 1 1; 1 0 1; 0 0 1; 0 0 -10];
costs = [2; 2; 3];

[energy, labels] = viterbiPottsMex(unary, costs);

verifyEqual(testCase, energy, -7);
verifyEqual(testCase, labels, [3; 3; 3; 3]);
end