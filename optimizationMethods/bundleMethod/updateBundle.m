function bundle = updateBundle(bundle, currentPoint, fValue, subgradient )
%updateBundle is an internal function for maximizeBundleMethodFixedSize
%
% Anton Osokin (firstname.lastname@gmail.com),  16.05.2013

% currentBundle = struct;
% currentBundle.maxBundleSize = options.maxBundleSize;
% currentBundle.numDimension = numVars;
% currentBundle.size = 0;
% currentBundle.points = zeros(currentBundle.numDimension, currentBundle.size);
% currentBundle.fValues = zeros(currentBundle.size, 1);
% currentBundle.subgradients = zeros(currentBundle.numDimension, currentBundle.size);
% currentBundle.pointSubgradProd = zeros(currentBundle.size, 1);
% currentBundle.grammMatrix = zeros(currentBundle.maxBundleSize, currentBundle.maxBundleSize);

    newBundlePos = bundle.size + 1;
    if newBundlePos > bundle.maxBundleSize % bundle will be too large
        % find a vector in the bundle to substitute
        
        curDot = double( currentPoint' * bundle.subgradients );
        bundleValues = bundle.fValues + curDot(:) - bundle.pointSubgradProd; 
        
        [~, newBundlePos] = max( bundleValues, [], 1 );
    end
    
    % add a vector to the bundle    
    bundle.size = max( bundle.size, newBundlePos);
%     bundle.points(:, newBundlePos) = currentPoint;
    bundle.fValues(newBundlePos) = fValue;
    bundle.subgradients(:, newBundlePos) = single( subgradient );
    bundle.pointSubgradProd(newBundlePos) = currentPoint' * subgradient;
    
    grammUpdate = double( subgradient' * bundle.subgradients );
    for iVector = 1 : bundle.size
        if iVector == newBundlePos
            bundle.grammMatrix( newBundlePos, newBundlePos ) = grammUpdate( newBundlePos );
        else
            bundle.grammMatrix( iVector, newBundlePos ) = grammUpdate( iVector );
            bundle.grammMatrix( newBundlePos, iVector ) = grammUpdate( iVector );
        end
    end
    
    
end
