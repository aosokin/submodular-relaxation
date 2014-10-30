function  gt = readMsrcGroundTruth( imageName )
% readMsrcGroundTruth is a helper function to work with the MSRC-21 dataset 
%
% See createDataset_MSRC for instructions
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

bmpImage = imread( imageName );

[ classColors, classNames ] = createMsrcClassColorTable;
% the first is void

gt = nan(size(bmpImage, 1), size(bmpImage, 2));

numClasses = length(classNames);
for iClass = 2 : numClasses
    curMask = true(size(bmpImage, 1), size(bmpImage, 2));
    for iChannel = 1 : 3
        curMask = curMask & (bmpImage(:, :, iChannel) == classColors(iClass, iChannel));
    end
    
    gt(curMask) = iClass - 1; % take case for the void class
end 	 	 	 	 

end

