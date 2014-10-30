function [ unary ] = readUnaryPotentials( fileName, imageHeight, imageWidth  )
% readUnaryPotentials is a helper function to work with the MSRC-21 dataset 
%
% See createDataset_MSRC for instructions
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

fid = fopen( fileName, 'r' );
if fid ~= 3
    error('readUnaryPotentials:badFile', ['Can not read ', fileName,]);
end

typeName = fscanf(fid, '%s', 1);
if ~strcmp(typeName, 'unary')
    fclose(fid);
    error('readUnaryPotentials:badFile', ['Error: cannot read the unary potentials from file ', fileName, ': the first word is wrong']);
end

numNodes = fscanf(fid, '%d', 1);
numLabels = fscanf(fid, '%d', 1);

[data, count] = fscanf(fid, '%f', numNodes * numLabels);
if count ~= numNodes * numLabels
    fclose(fid);    
    error('readUnaryPotentials:badFile', ['Error: cannot read the unary potentials from file ', fileName, ': wrong number of numbers']);
end
fclose(fid);

if numNodes ~= imageHeight * imageWidth
    error('readUnaryPotentials:badFile', 'Error: image size is not compatible with the unary potentials')
end
    
unary = reshape(data, [numLabels, numNodes] );

% reorder the nodes to the column-first ordering
imageIds_correct = reshape(1 : imageHeight * imageWidth, [imageHeight, imageWidth] );
imageIds_current = imageIds_correct(end:-1:1, :);
imageIds_current = imageIds_current';
[~, inverseIds] = sort(imageIds_current(:));
unary = unary(:, inverseIds);

end

