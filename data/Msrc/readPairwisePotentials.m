function [ pairwise ] = readPairwisePotentials( fileName, imageHeight, imageWidth  )
% readPairwisePotentials is a helper function to work with the MSRC-21 dataset 
%
% See createDataset_MSRC for instructions
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

fid = fopen( fileName, 'r' );
if fid ~= 3
    error('readPairwisePotentials:badFile', ['Can not read ', fileName,]);
end

typeName = fscanf(fid, '%s', 1);
if ~strcmp(typeName, 'pairwisePotts')
    fclose(fid);
    error('readPairwisePotentials:badFile', ['Error: cannot read the pairwise potentials from file ', fileName, ': the first word is wrong']);
end

numEdges = fscanf(fid, '%d', 1);

[data, count] = fscanf(fid, '%f', numEdges * 3);
if count ~= numEdges * 3
    fclose(fid);    
    error('readPairwisePotentials:badFile', ['Error: cannot read the pairwise potentials from file ', fileName, ': wrong number of numbers']);
end
fclose(fid);

pairwise = reshape( data, [3, numEdges] )';

% reorder the nodes to the column-first ordering
imageIds_correct = reshape(1 : imageHeight * imageWidth, [imageHeight, imageWidth] );
imageIds_current = imageIds_correct(end : -1 : 1, :);
imageIds_current = imageIds_current';
newIds = imageIds_current(:);

pairwise(:, 1 : 2) = newIds( pairwise(:, 1 : 2) + 1 );

toSwap = pairwise(:, 1) > pairwise(:, 2);

tmp = pairwise(toSwap, 1);
pairwise(toSwap, 1) = pairwise(toSwap, 2);
pairwise(toSwap, 2) = tmp;

pairwise = sparse( pairwise(:, 1), pairwise(:, 2), pairwise(:, 3), imageHeight * imageWidth, imageHeight * imageWidth ); 

end

