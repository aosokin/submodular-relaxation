function hammingError = computeHammingError( groundTruth, result )
% rcomputeHammingError is a helper function to work with the MSRC-21 dataset 
%
% See createDataset_MSRC for instructions
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

if numel(size(groundTruth))~= 2 || numel(size(groundTruth))~= 2 || any( size(groundTruth) ~= size(result) )
    error('computeHammingError:wrongInput', 'Error: the ground truth and the result are not compatible');
end

groundTruth = groundTruth(:);
result = result(:);

goodNodes = ~isnan(groundTruth);

groundTruth = groundTruth( goodNodes );
result = result( goodNodes );

hammingError = sum(groundTruth ~= result) / length(groundTruth);

end

