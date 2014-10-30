function experiment_MSRC(datasetId, matlabDatasetFolder, imageIds, matlabResultFolder)
% experiment_MSRC conducts the experiment section 7.2 of the journal paper, results: fig. 5
%
% datasetId: can be '1seg' or '3seg'; (default: '1seg')
% matlabDatasetFolder: path to the dataset in MATLAB format (see createDataset_MSRC for details) (default: data/Msrc/MatlabDataset)
% imageIds: indices of images for which to compute the results (default: all images)
% matlabResultFolder: folder to store results (default: data/Msrc/MatlabResults_{datasetId})
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

setup_SMR;

%% parse the parameters
if ~exist('datasetId', 'var') || isempty(datasetId)
    datasetId = '1seg';
end
if ~isequal(datasetId, '1seg') && ~isequal(datasetId, '3seg')
    error('experiment_MSRC:badDatasetId', 'Bad value of datasetId, should be 1seg or 3seg');
end

if ~exist('matlabDatasetFolder', 'var') || isempty(matlabDatasetFolder)
    matlabDatasetFolder = 'data/Msrc/MatlabDataset';
end
if ~ischar(matlabDatasetFolder) || ~exist(matlabDatasetFolder, 'dir') 
    error('experiment_MSRC:badMatlabDatasetFolder', 'matlabDatasetFolder should contain the path to the data (output of createDataset_MSRC)');
end

if ~exist('matlabResultFolder', 'var') || isempty(matlabResultFolder)
    matlabResultFolder = 'data/Msrc';
end
if ~ischar(matlabResultFolder) 
    error('experiment_MSRC:badMatlabResultsFolder', 'matlabResultsFolder should be a string');
end
matlabResultFolder = fullfile(matlabResultFolder, ['MatlabResults_', datasetId]);
if ~exist(matlabResultFolder, 'dir')
    mkdir(matlabResultFolder);
end

%% get the list of images
listOfFiles = dir(fullfile(matlabDatasetFolder, '*.mat'));
numImages = length(listOfFiles);

if ~exist('imageIds', 'var') || isempty(imageIds)
    imageIds = 1 : numImages;
end

%% read all image data
resultImageNames = cell(0,0);
for iImage = intersect(1 : numImages, imageIds);
    [~, curImageName, ~] = fileparts(listOfFiles(iImage).name);
    resultImageNames{end + 1} = curImageName;
    fprintf( 'Working with image %d of %d: %s\n', iImage, numImages, curImageName );
    
    s = RandStream('mt19937ar','Seed',1);
    RandStream.setGlobalStream(s);
    
    load( fullfile(matlabDatasetFolder, [curImageName, '.mat']), 'imageName', 'unary', 'pairwisePotts', 'rgbImage', 'groundTruth', ...
                ['highOrderNodes_', datasetId], ['highOrderParameters_', datasetId] );
    eval(['highOrderNodes = ', 'highOrderNodes_', datasetId, ';']);
    eval(['highOrderParameters = ', 'highOrderParameters_', datasetId, ';']);
            
    if ~isequal(imageName, curImageName)
        error('experiment_MSRC:badImageName', ['Image file ', curImageName, ' is not consistent!' ]);
    end

    unary = double(unary);
    imageWidth = size(groundTruth, 2);
    imageHeight = size(groundTruth, 1);
    numNodes = size(unary, 2);
    numLabels = size(unary, 1);
    
    [~, labels_unaryOnly] = min(unary, [], 1);
    energy_unaryOnly = computeEnergy_highOrderPotts(unary, pairwisePotts, labels_unaryOnly(:), highOrderNodes, highOrderParameters); 
    hammingError_unaryOnly = computeHammingError(groundTruth, reshape(labels_unaryOnly, [imageHeight, imageWidth]) );
    
    [labels_unaryPairwise] = alphaExpansionRobustHighOrderPottsMex(unary, pairwisePotts);
    energy_unaryPairwise = computeEnergy_highOrderPotts(unary, pairwisePotts, labels_unaryPairwise(:), highOrderNodes, highOrderParameters);
    hammingError_unaryPairwise = computeHammingError(groundTruth, reshape(labels_unaryPairwise, [imageHeight, imageWidth]) );
    
    [labels_unaryPairwiseHo, energy, energyPlot, timePlot] = alphaExpansionRobustHighOrderPottsMex(unary, pairwisePotts, highOrderNodes, highOrderParameters);
    energy_unaryPairwiseHo = computeEnergy_highOrderPotts(unary, pairwisePotts, labels_unaryPairwiseHo(:), highOrderNodes, highOrderParameters);
    if abs(energy - energy_unaryPairwiseHo) > 1e-5
        warning('Something is wrong with the energy computation: energy returned by alphaExpansionRobustHighOrderPottsMex is not equal to energy computed by computeEnergy_highOrderPotts\n');
    end
    hammingError_unaryPairwiseHo = computeHammingError(groundTruth, reshape(labels_unaryPairwiseHo, [imageHeight, imageWidth]) );
    timeAlphaExp = max(timePlot);
    
    % run SMD
    computeSmrDualDynamic_highOrderPotts_clearGlobal();
    options = struct;
    options.maxIter = 300;
    options.argTol = 1e-2;
    options.funcGetPrimalLabeling = @(curLabels) deal(computeEnergy_highOrderPotts(unary, pairwisePotts, curLabels, highOrderNodes, highOrderParameters), curLabels);
    oracle = @(dualVars) computeSmrDualDynamic_highOrderPotts(unary, pairwisePotts, dualVars, highOrderNodes, highOrderParameters);

    [ bestDualPoint_smr, bestDual_smr, bestPrimal_smr, bestLabeling_smr, timePlot_smr, dualPlot_smr, primalPlot_smr ] = maximizeHanso( numNodes, oracle, options );
    hammingError_unaryPairwiseHo_smr = computeHammingError(groundTruth, reshape(bestLabeling_smr, [imageHeight, imageWidth]) );
    computeSmrDualDynamic_highOrderPotts_clearGlobal();
    timeSmr = max(timePlot_smr);   
    
    % run CWD
    options = struct;
    options.maxIter = 300;
    options.argTol = 1e-2;
    options.funcGetPrimalLabeling = @(curLabels) deal(computeEnergy_highOrderPotts(unary, pairwisePotts, curLabels, highOrderNodes, highOrderParameters), curLabels);
    [pairwiseDirectionalCosts] = separatePairwiseToDirections( pairwisePotts, [imageHeight; imageWidth] );
    oracle = @(dualVars) computeDdtrwDual_highOrderPotts(unary, pairwiseDirectionalCosts, dualVars, highOrderNodes, highOrderParameters);
    
    numDualVars = numNodes * numLabels;
    if isequal(datasetId, '1seg') 
        numDualVars = numNodes * numLabels * 4;
    end
    if isequal(datasetId, '3seg')
        numDualVars = numNodes * numLabels * 6;
    end
    
    [ bestDualPoint_cwd, bestDual_cwd, bestPrimal_cwd, bestLabeling_cwd, timePlot_cwd, dualPlot_cwd, primalPlot_cwd ] = maximizeHanso( numDualVars, oracle, options );
    hammingError_unaryPairwiseHo_cwd = computeHammingError(groundTruth, reshape(bestLabeling_cwd, [imageHeight, imageWidth]) );
    timeCwd = max(timePlot_cwd);   
    
    save( fullfile( matlabResultFolder, [curImageName, '.mat']), 'imageWidth', 'imageHeight', 'timeAlphaExp', 'timeSmr', 'timeCwd', ...   
     'matlabDatasetFolder', 'labels_unaryOnly', 'energy_unaryOnly',  'hammingError_unaryOnly', ...
     'labels_unaryPairwise', 'energy_unaryPairwise', 'hammingError_unaryPairwise', ...
     'labels_unaryPairwiseHo', 'energy_unaryPairwiseHo', 'hammingError_unaryPairwiseHo',...
     'bestDualPoint_smr', 'bestDual_smr', 'bestPrimal_smr', 'bestLabeling_smr', 'timePlot_smr', 'dualPlot_smr', 'primalPlot_smr', 'hammingError_unaryPairwiseHo_smr', ...
     'bestDualPoint_cwd', 'bestDual_cwd', 'bestPrimal_cwd', 'bestLabeling_cwd', 'timePlot_cwd', 'dualPlot_cwd', 'primalPlot_cwd', 'hammingError_unaryPairwiseHo_cwd');
end

makePlots_MSRC(datasetId, matlabResultFolder, resultImageNames);

end
