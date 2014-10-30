function matlabDatasetFolder = createDataset_MSRC
% createDataset_MSRC reads the data from the MSRC dataset and constructs the features used in the journal paper.
%
% The MSRC dataset: http://www.inf.ethz.ch/personal/ladickyl/Msrc.zip
% The original source of the ALE library: http://www.inf.ethz.ch/personal/ladickyl/ALE.zip
%
% If you are using this code and data please cite the following papers:
%
%   1) Pushmeet Kohli, Lubor Ladicky, Philip H. S. Torr.
%   Robust Higher Order Potentials for Enforcing Label Consistency
%   International Journal of Computer Vision, 82(3):302-324, 2009.
%
%   2) Jamie Shotton, John Winn, Carsten Rother, Antonio Criminisi
%   TextonBoost: Joint Appearance, Shape and Context Modeling for Multi-class Object Recognition and Segmentation
%   Computer Vision – ECCV 2006
%   Lecture Notes in Computer Science Volume 3951, 2006, pp 1-15
%
% createDataset_MSRC relies on the ALE library (for features) which works under Windows only.
%
% Installation instructions:
% 1) Unpack the ./data/ALE_changed.zip to ./data (Note that this code was written by Lubor Ladicky and is released under research0only license).
%       This version was slightly modified by Anton Osokin to explicitly save the potentials.
% 2) Compile the code (tested with MSVC 2010 ans MSVC 2012). ./data/ALE_binaries.zip contains the Windows x64 binaries compiled with MSVC 2012 that might work for you.
% 3) Put the ALE binary file together with the required dlls (ILU.dll, ILUT.dll, DevIL.dll) to ./data/
% 4) Download the MSRC-21 dataset from http://www.inf.ethz.ch/personal/ladickyl/Msrc.zip and put it to ./data/Msrc folder in such a way that folder ./data/Msrc/Images contains images of the dataset
% 5) Check that files ./data/Msrc/Test.txt and ./data/Msrc/Train.txt exist. If not dowload them from http://jamie.shotton.org/work/data/TextonBoostSplits.zip
% 6) After all this setup run ./data/ALE.exe to compute the features. Note that this will take some time.
% 7) Run createDataset_MSRC to create dataset suitable for submodular-relaxation framework
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

%% initialization
%msrcDatasetFolder = 'Msrc';
msrcDatasetFolder = '/local/aosokin/data/Msrc';

potentialsFolder = fullfile(msrcDatasetFolder, 'Result', 'Potentials');
imagesFolder = fullfile(msrcDatasetFolder, 'Images');
groundtruthFolder = fullfile(msrcDatasetFolder, 'GroundTruth');

segmentationFolder_3seg = cell(3, 1);
segmentationFolder_3seg{1} = fullfile(msrcDatasetFolder, 'Result/MeanShift/70x65');
segmentationFolder_3seg{2} = fullfile(msrcDatasetFolder, 'Result/MeanShift/70x95');
segmentationFolder_3seg{3} = fullfile(msrcDatasetFolder, 'Result/MeanShift/70x145');

segmentationFolder_1seg = cell(1, 1);
segmentationFolder_1seg{1} = fullfile(msrcDatasetFolder, 'Result/MeanShift/70x95');


matlabDatasetFolder = fullfile( msrcDatasetFolder, 'MatlabDataset');
if ~exist(matlabDatasetFolder, 'dir')
    mkdir(matlabDatasetFolder);
end

% parameters of the high-order potentials
options.theta_alpha = 0.8;
options.theta_hp = 0.2;
options.theta_hv = 0.5;
options.theta_hb = 12;
options.Q_ratio = 0.1;

%% get the list of images
imageListFile = fullfile(msrcDatasetFolder, 'Test.txt');
fileID = fopen(imageListFile, 'r');
if fileID ~= 3
    error('createDataset_MSRC:testTxtMissing', ['File ', imageListFile,' could not be opened']);
end
lines = textscan(fileID, '%s');
fclose(fileID);
lines = lines{1};
numImages = length(lines);
imageNames = cell( numImages, 1 );
for iImage = 1 : numImages
    [~, fileName, ~] = fileparts( lines{iImage} );
    imageNames{ iImage } = fileName;
end
clear lines fileID fileName

%% Prepare the MATLAB dataset
for iImage = 1 : length(imageNames)
    imageName = imageNames{iImage};
    
    fprintf('Working with image %d of %d: %s\n', iImage, length(imageNames), imageName);
    
    rgbImage = imread( fullfile(imagesFolder, [imageName,'.bmp']) ) ;
    groundTruth = readMsrcGroundTruth( fullfile(groundtruthFolder, [imageName,'.bmp']) );
    
    segments_3seg = cell(length(segmentationFolder_3seg), 1);
    for iSegmentation = 1 : length(segmentationFolder_3seg)
        segments_3seg{iSegmentation} = readBinaryFileAle( fullfile(segmentationFolder_3seg{iSegmentation}, [imageName, '.msh']) );
    end
    
    segments_1seg = cell(length(segmentationFolder_1seg), 1);
    for iSegmentation = 1 : length(segmentationFolder_1seg)
        segments_1seg{iSegmentation} = readBinaryFileAle( fullfile(segmentationFolder_1seg{iSegmentation}, [imageName, '.msh']) );
    end
    
    
    imageWidth = size(rgbImage, 2);
    imageHeight = size(rgbImage, 1);
    unary = readUnaryPotentials( fullfile(potentialsFolder, [imageName,'_0.data']), imageHeight, imageWidth );
    unary = single(unary);
    pairwisePotts = readPairwisePotentials( fullfile(potentialsFolder, [imageName,'_1.data']), imageHeight, imageWidth );
    [highOrderNodes_3seg, highOrderParameters_3seg ] = createHighOrderPotentials( rgbImage, segments_3seg, options );
    [highOrderNodes_1seg, highOrderParameters_1seg ] = createHighOrderPotentials( rgbImage, segments_1seg, options );
    
    save( fullfile(matlabDatasetFolder, [imageName, '.mat']), 'imageName', 'unary', 'pairwisePotts', 'rgbImage', 'groundTruth', ...
                'highOrderNodes_3seg', 'highOrderParameters_3seg', 'highOrderNodes_1seg', 'highOrderParameters_1seg');
end

end









































