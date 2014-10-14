function createDataset_NSMR
% createDataset_NSMR creates artificial dataset for the NSMR experiment of the journal paper
%
% Anton Osokin (firstname.lastname@gmail.com),  14.10.2014

%% dataset folder
curFile = mfilename('fullpath');
dataPath = fileparts(curFile);

%% dataset parameters
nObjects = 50;
nRows = 20;
nCols = 20;
nLabels = 5;
weightPairwise = sqrt(0.5);

%% fix rand seed
s = RandStream('mt19937ar', 'Seed', 1);
RandStream.setGlobalStream(s);

%% create data
dataset = cell(nObjects, 1);
for iObject = 1 : nObjects
    dataset{iObject}.type = 'toy';
    dataset{iObject}.name = ['nsmr', num2str(iObject, '%03d')];
    dataset{iObject}.nRows = nRows;
    dataset{iObject}.nCols = nCols;
    dataset{iObject}.nLabels = nLabels;
    
    dataset{iObject}.unary = randn(nLabels, nRows * nCols);
    pottsWeight = struct;
    
    pottsWeight.vertCost = randn(nRows - 1, nCols) * weightPairwise;
    pottsWeight.horCost = randn(nRows, nCols - 1) * weightPairwise;
    dataset{iObject}.pairwisePotts = buildNeighborhoodGridPotts(nRows, nCols, pottsWeight, '4');
end

%% save dataset
save(fullfile(dataPath, 'dataset_NSMR.mat'), 'dataset');

end
