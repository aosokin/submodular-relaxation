function data = readBinaryFileAle( fileName, dataType )
% readBinaryFileAle is a helper function to work with the MSRC-21 dataset 
%
% See createDataset_MSRC for instructions
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

if ~exist('dataType', 'var')
    dataType = 'int32';
end

fid = fopen( fileName, 'r');

if fid ~= 3
    error('readBinaryFileAle:badFile', ['Can not read ', fileName,]);
end

width = fread(fid, 1, 'int32');
height = fread(fid, 1, 'int32');
numChannels = fread(fid, 1, 'int32');
targetCount = width * height * numChannels;
[data, countRead] = fread(fid, targetCount, dataType);

fclose(fid);

if countRead ~= targetCount
    error(['Error while reading file ', fileName]);
end
data = double(data);
data = reshape(data, [numChannels, width, height]);
data = permute(data, [3 2 1]);
data = data(end:-1:1, :, :);

end

