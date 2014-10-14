function createDataset_alahariPAMI10
% createDataset_alahariPAMI10 reads data provided by K. Alahari and transforms 
% all instances to the internal format. 
%
% Link to get the data: http://www.di.ens.fr/~alahari/data/pami10data.tgz
%
% If you are using the data please cite the following paper in any resulting publication:
% K. Alahari, P. Kohli, and P. H. S. Torr, 
% Dynamic hybrid algorithms for MAP inference in discreteMRFs, 
% IEEE TPAMI, vol. 32, no. 10, pp. 1846-1857, 2010.
%
% Anton Osokin (firstname.lastname@gmail.com),  14.10.2014

dataUrl = 'http://www.di.ens.fr/~alahari/data/pami10data.tgz';

% dataset folder
curFile = mfilename('fullpath');
dataPath = fileparts(curFile);
dataFileName = 'pami10data.tgz';
if ~exist( fullfile(dataPath, dataFileName) , 'file' )
    fprintf('Downloading dataset from %s\n', dataUrl);
    urlwrite( dataUrl, fullfile('data', dataFileName) );
end
fprintf('Unpacking the dataset\n');
untar(fullfile(dataPath, dataFileName), dataPath);

dataDir = fullfile(dataPath, 'pami10data');
dataset = cell(0, 0);

fprintf('Creating instaces\n');
%% stereo data
stereoDir = fullfile(dataDir, 'stereo');

% tsukuba
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'stereo';
dataset{iEnergy}.name = 'tsukuba';

dataset{iEnergy}.nRows = 288;
dataset{iEnergy}.nCols = 384;
dataset{iEnergy}.nLabels = 16;

dataset{iEnergy}.pottsCoefficient = 20;

% read and reformat unary potentials
curUnary = dlmread(fullfile(stereoDir, 'tsukuba_datacost.txt'));
curUnary = reshape(curUnary, [dataset{iEnergy}.nCols, dataset{iEnergy}.nRows, dataset{iEnergy}.nLabels]);
curUnary = permute(curUnary, [2, 1, 3]); % change dimensions to have standard appearance
dataset{iEnergy}.unary = reshape( curUnary, [], dataset{iEnergy}.nLabels )';
dataset{iEnergy}.pairwisePotts = buildNeighborhoodGridPotts(dataset{iEnergy}.nRows, dataset{iEnergy}.nCols, dataset{iEnergy}.pottsCoefficient, '4');

% venus
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'stereo';
dataset{iEnergy}.name = 'venus';

dataset{iEnergy}.nRows = 383;
dataset{iEnergy}.nCols = 434;
dataset{iEnergy}.nLabels = 20;

dataset{iEnergy}.pottsCoefficient = 20;

% read and reformat unary potentials
curUnary = dlmread(fullfile(stereoDir, 'venus_datacost.txt'));
curUnary = reshape(curUnary, [dataset{iEnergy}.nCols, dataset{iEnergy}.nRows, dataset{iEnergy}.nLabels]);
curUnary = permute(curUnary, [2, 1, 3]); % change dimensions to have standard appearance
dataset{iEnergy}.unary = reshape( curUnary, [], dataset{iEnergy}.nLabels )';
dataset{iEnergy}.pairwisePotts = buildNeighborhoodGridPotts(dataset{iEnergy}.nRows, dataset{iEnergy}.nCols, dataset{iEnergy}.pottsCoefficient, '4');

% cones
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'stereo';
dataset{iEnergy}.name = 'cones';

dataset{iEnergy}.nRows = 375;
dataset{iEnergy}.nCols = 450;
dataset{iEnergy}.nLabels = 60;

dataset{iEnergy}.pottsCoefficient = 20;

% read and reformat unary potentials
curUnary = dlmread(fullfile(stereoDir, 'cones_datacost.txt'));
curUnary = reshape(curUnary, [dataset{iEnergy}.nCols, dataset{iEnergy}.nRows, dataset{iEnergy}.nLabels]);
curUnary = permute(curUnary, [2, 1, 3]); % change dimensions to have standard appearance
dataset{iEnergy}.unary = reshape( curUnary, [], dataset{iEnergy}.nLabels )';
dataset{iEnergy}.pairwisePotts = buildNeighborhoodGridPotts(dataset{iEnergy}.nRows, dataset{iEnergy}.nCols, dataset{iEnergy}.pottsCoefficient, '4');

% teddy
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'stereo';
dataset{iEnergy}.name = 'teddy';

dataset{iEnergy}.nRows = 375;
dataset{iEnergy}.nCols = 450;
dataset{iEnergy}.nLabels = 60;

dataset{iEnergy}.pottsCoefficient = 20;

% read and reformat unary potentials
curUnary = dlmread(fullfile(stereoDir, 'teddy_datacost.txt'));
curUnary = reshape(curUnary, [dataset{iEnergy}.nCols, dataset{iEnergy}.nRows, dataset{iEnergy}.nLabels]);
curUnary = permute(curUnary, [2, 1, 3]); % change dimensions to have standard appearance
dataset{iEnergy}.unary = reshape( curUnary, [], dataset{iEnergy}.nLabels )';
dataset{iEnergy}.pairwisePotts = buildNeighborhoodGridPotts(dataset{iEnergy}.nRows, dataset{iEnergy}.nCols, dataset{iEnergy}.pottsCoefficient, '4');

%% Object segmentation
objSegDir = fullfile(dataDir, 'object-seg');

% 349
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'obj-seg';
dataset{iEnergy}.name = '349';

dataset{iEnergy}.nRows = 213;
dataset{iEnergy}.nCols = 320;
dataset{iEnergy}.nLabels = 4;

dataFile = dlmread(fullfile(objSegDir, 'energy-349.txt'));
curUnary = dataFile( 2 : dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 1, 2 : 1 + dataset{iEnergy}.nLabels);
curNeighbors = dataFile( dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 2 : end, 1 : 3);

dataset{iEnergy}.unary = curUnary';
dataset{iEnergy}.pairwisePotts = sparse( [ curNeighbors(:, 1) + 1; curNeighbors(:, 2) + 1], [curNeighbors(:, 2) + 1; curNeighbors(:, 1) + 1], [curNeighbors(:, 3); curNeighbors(:, 3)] ); 

% 353
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'obj-seg';
dataset{iEnergy}.name = '353';

dataset{iEnergy}.nRows = 213;
dataset{iEnergy}.nCols = 320;
dataset{iEnergy}.nLabels = 7;

dataFile = dlmread(fullfile(objSegDir, 'energy-353.txt'));
curUnary = dataFile( 2 : dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 1, 2 : 1 + dataset{iEnergy}.nLabels);
curNeighbors = dataFile( dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 2 : end, 1 : 3);

dataset{iEnergy}.unary = curUnary';
dataset{iEnergy}.pairwisePotts = sparse( [ curNeighbors(:, 1) + 1; curNeighbors(:, 2) + 1], [curNeighbors(:, 2) + 1; curNeighbors(:, 1) + 1], [curNeighbors(:, 3); curNeighbors(:, 3)] ); 

% 353 8-connected
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'obj-seg';
dataset{iEnergy}.name = '353-8con';

dataset{iEnergy}.nRows = 320;
dataset{iEnergy}.nCols = 213;
dataset{iEnergy}.nLabels = 7;

dataFile = dlmread(fullfile(objSegDir, 'energy-353-8neighbourhood.txt'));
curUnary = dataFile( 2 : dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 1, 2 : 1 + dataset{iEnergy}.nLabels);
curNeighbors = dataFile( dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 2 : end, 1 : 3);

dataset{iEnergy}.unary = curUnary';
dataset{iEnergy}.pairwisePotts = sparse( [ curNeighbors(:, 1) + 1; curNeighbors(:, 2) + 1], [curNeighbors(:, 2) + 1; curNeighbors(:, 1) + 1], [curNeighbors(:, 3); curNeighbors(:, 3)] ); 

% 358
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'obj-seg';
dataset{iEnergy}.name = '358';

dataset{iEnergy}.nRows = 213;
dataset{iEnergy}.nCols = 320;
dataset{iEnergy}.nLabels = 5;

dataFile = dlmread(fullfile(objSegDir, 'energy-358.txt'));
curUnary = dataFile( 2 : dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 1, 2 : 1 + dataset{iEnergy}.nLabels);
curNeighbors = dataFile( dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 2 : end, 1 : 3);

dataset{iEnergy}.unary = curUnary';
dataset{iEnergy}.pairwisePotts = sparse( [ curNeighbors(:, 1) + 1; curNeighbors(:, 2) + 1], [curNeighbors(:, 2) + 1; curNeighbors(:, 1) + 1], [curNeighbors(:, 3); curNeighbors(:, 3)] ); 

% 416
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'obj-seg';
dataset{iEnergy}.name = '416';

dataset{iEnergy}.nRows = 213;
dataset{iEnergy}.nCols = 320;
dataset{iEnergy}.nLabels = 5;

dataFile = dlmread(fullfile(objSegDir, 'energy-416.txt'));
curUnary = dataFile( 2 : dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 1, 2 : 1 + dataset{iEnergy}.nLabels);
curNeighbors = dataFile( dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 2 : end, 1 : 3);

dataset{iEnergy}.unary = curUnary';
dataset{iEnergy}.pairwisePotts = sparse( [ curNeighbors(:, 1) + 1; curNeighbors(:, 2) + 1], [curNeighbors(:, 2) + 1; curNeighbors(:, 1) + 1], [curNeighbors(:, 3); curNeighbors(:, 3)] ); 

% 552
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'obj-seg';
dataset{iEnergy}.name = '552';

dataset{iEnergy}.nRows = 213;
dataset{iEnergy}.nCols = 320;
dataset{iEnergy}.nLabels = 8;

dataFile = dlmread(fullfile(objSegDir, 'energy-552.txt'));
curUnary = dataFile( 2 : dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 1, 2 : 1 + dataset{iEnergy}.nLabels);
curNeighbors = dataFile( dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 2 : end, 1 : 3);

dataset{iEnergy}.unary = curUnary';
dataset{iEnergy}.pairwisePotts = sparse( [ curNeighbors(:, 1) + 1; curNeighbors(:, 2) + 1], [curNeighbors(:, 2) + 1; curNeighbors(:, 1) + 1], [curNeighbors(:, 3); curNeighbors(:, 3)] ); 
 

%% Colour segmentation
objSegDir = fullfile(dataDir, 'colour-seg');

% garden
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'colour-seg';
dataset{iEnergy}.name = 'garden';

dataset{iEnergy}.nRows = 175;
dataset{iEnergy}.nCols = 120;
dataset{iEnergy}.nLabels = 4;

dataFile = dlmread(fullfile(objSegDir, 'garden4_costs.txt'));
curUnary = dataFile( 2 : dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 1, 2 : 1 + dataset{iEnergy}.nLabels);
curNeighbors = dataFile( dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 2 : end, 1 : 3);

dataset{iEnergy}.unary = curUnary';
dataset{iEnergy}.pairwisePotts = sparse( [ curNeighbors(:, 1) + 1; curNeighbors(:, 2) + 1], [curNeighbors(:, 2) + 1; curNeighbors(:, 1) + 1], [curNeighbors(:, 3); curNeighbors(:, 3)] ); 

% cow3
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'colour-seg';
dataset{iEnergy}.name = 'cow3';

dataset{iEnergy}.nRows = 720;
dataset{iEnergy}.nCols = 576;
dataset{iEnergy}.nLabels = 3;

dataFile = dlmread(fullfile(objSegDir, 'cow3_costs.txt'));
curUnary = dataFile( 2 : dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 1, 2 : 1 + dataset{iEnergy}.nLabels);
curNeighbors = dataFile( dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 2 : end, 1 : 3);

dataset{iEnergy}.unary = curUnary';
dataset{iEnergy}.pairwisePotts = sparse( [ curNeighbors(:, 1) + 1; curNeighbors(:, 2) + 1], [curNeighbors(:, 2) + 1; curNeighbors(:, 1) + 1], [curNeighbors(:, 3); curNeighbors(:, 3)] ); 

% cow4
iEnergy = length(dataset) + 1;
dataset{iEnergy} = struct;
dataset{iEnergy}.type = 'colour-seg';
dataset{iEnergy}.name = 'cow4';

dataset{iEnergy}.nRows = 720;
dataset{iEnergy}.nCols = 576;
dataset{iEnergy}.nLabels = 4;

dataFile = dlmread(fullfile(objSegDir, 'cow4_costs.txt'));
curUnary = dataFile( 2 : dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 1, 2 : 1 + dataset{iEnergy}.nLabels);
curNeighbors = dataFile( dataset{iEnergy}.nRows * dataset{iEnergy}.nCols + 2 : end, 1 : 3);

dataset{iEnergy}.unary = curUnary';
dataset{iEnergy}.pairwisePotts = sparse( [ curNeighbors(:, 1) + 1; curNeighbors(:, 2) + 1], [curNeighbors(:, 2) + 1; curNeighbors(:, 1) + 1], [curNeighbors(:, 3); curNeighbors(:, 3)] ); 

dataset = dataset(:);
 
%% normalize instances

problemSet =  1 : length(dataset);
alphaExpRenormalizationValue = 100;
roundingPrecision = 1e-3; 
    
for iEnergy = problemSet
    fprintf('Normalizing energy %d\n', iEnergy);
    
    %% round potentials
    curUnary = round(dataset{iEnergy}.unary / roundingPrecision) * roundingPrecision;
    curPairwise = round(dataset{iEnergy}.pairwisePotts / roundingPrecision) * roundingPrecision;
    
    %% run alpha-expansion
    [alphaExp_labels, alphaExp_energy, alphaExp_time] = alphaExpansionPotts(curUnary, curPairwise, roundingPrecision);
    
    alphaExp_trueEnergy = computeEnergyPotts(curUnary, curPairwise, alphaExp_labels);
    if abs(alphaExp_trueEnergy - alphaExp_energy) > 1e-6 * (abs(max(alphaExp_trueEnergy, alphaExp_energy)) + 1e-2)
        warning(['Problem: ', num2str(iEnergy), ' alpha-exp energy output does not equal recomputed energy, difference: ', num2str(abs(alphaExp_trueEnergy - alphaExp_energy)) ]);
    end
    
    %% normalize
	bestUnaryValue = min(curUnary, [], 1);
    
    alphaExp_energy_new = alphaExp_energy - sum(bestUnaryValue, 2);
    
    newUnary = bsxfun(@minus, curUnary, bestUnaryValue) / alphaExp_energy_new * alphaExpRenormalizationValue;
    newPairwise = curPairwise / alphaExp_energy_new * alphaExpRenormalizationValue;
    
    normalizationMinus = sum(bestUnaryValue, 2);
    normalizationDivide = alphaExp_energy_new / alphaExpRenormalizationValue;

    dataset{iEnergy}.unary = newUnary ;
    dataset{iEnergy}.pairwisePotts = newPairwise;
    dataset{iEnergy}.normalizationMinus = normalizationMinus;
    dataset{iEnergy}.normalizationDivide = normalizationDivide;
end

%% saving results
fprintf('Save dataset to %s\n', fullfile(dataPath, 'dataset_alahariPAMI10.mat'));
save(fullfile(dataPath, 'dataset_alahariPAMI10.mat'), 'dataset');

%% clean up
rmdir(dataDir, 's')

end




















