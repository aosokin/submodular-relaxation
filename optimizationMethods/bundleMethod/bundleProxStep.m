function [ nextPoint, bundleValue] = bundleProxStep( bundle, currentMeanPoint, inverseStepSize)
%bundleProxStep is an internal function for maximizeBundleMethodFixedSize
%
% Anton Osokin (firstname.lastname@gmail.com),  16.05.2013

% bundle = struct;
% bundle.maxBundleSize = options.maxBundleSize;
% bundle.numDimension = numVars;
% bundle.size = 0;
% bundle.points = zeros(currentBundle.numDimension, currentBundle.size);
% bundle.fValues = zeros(currentBundle.size, 1);
% bundle.subgradients = zeros(currentBundle.numDimension, currentBundle.size);
% bundle.pointSubgradProd = zeros(currentBundle.size, 1);
% bundle.grammMatrix = zeros(currentBundle.maxBundleSize, currentBundle.maxBundleSize);


% construct a QP problem

bundleFValues = bundle.fValues( 1 : bundle.size );
bundlePointSubgradProd = bundle.pointSubgradProd( 1 : bundle.size );

curDot = currentMeanPoint' * bundle.subgradients;
curDot = double(curDot( 1 : bundle.size ));
f = bundleFValues + curDot(:) - bundlePointSubgradProd; 

f = f';
H = (bundle.grammMatrix) / (inverseStepSize);
H = H(1 : bundle.size, 1 : bundle.size);
Aeq = ones( 1, bundle.size);
beq = 1;
lb = zeros( bundle.size, 1 );
ub = ones( bundle.size, 1 );
options = optimset( 'Algorithm', 'interior-point-convex', 'Display', 'off' );

[xi, fValue, exitFlag] = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options);

if exitFlag ~= 1
    save('bundleOptimizationError.mat');
    warning('bundleProxStep:quadprog', ['optimization exited with exit code ', num2str(exitFlag)]);
end

xiTemp = [xi; zeros(bundle.maxBundleSize - bundle.size, 1)] / inverseStepSize;
curProd = double( bundle.subgradients * xiTemp );
nextPoint = curProd + currentMeanPoint;

curDot = nextPoint' * bundle.subgradients;
curDot = double( curDot( 1 : bundle.size ) );
bundleValue = min( bundleFValues + curDot(:) - bundlePointSubgradProd, [], 1);

% check duality
computedFValue = bundleValue - inverseStepSize / 2 * norm(nextPoint - currentMeanPoint)^2;
if abs( fValue - computedFValue) > 1e-2
%     save('bundleOptimizationError.mat');
    warning('bundleProxStep:duality', ['Primal and dual values do not coinside, diff: ', num2str(abs( fValue - computedFValue))]);
end
    
end
