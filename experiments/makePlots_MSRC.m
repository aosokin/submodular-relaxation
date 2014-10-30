function makePlots_MSRC(datasetId, matlabResultFolder, resultImageNames)
% experiment_MSRC conducts the experiment section 7.2 of the journal paper, results: fig. 5, (a) - 1seg; (b) - 3seg
%
% datasetId: can be '1seg' or '3seg'; (default: '1seg')
% matlabResultFolder: folder to store results (default: data/Msrc/MatlabResults_{datasetId})
% resultImageNames: cell array of image names for which to compute the results, file fullfile(matlabResultFolder, [resultsImageNames{i},'.mat']) has to be present
%       (default: all the files in a directory)
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

setup_SMR;

%% parse the parameters
if ~exist('datasetId', 'var') || isempty(datasetId)
    datasetId = '1seg';
end
if ~isequal(datasetId, '1seg') && ~isequal(datasetId, '3seg')
    error('makePlots_MSRC:badDatasetId', 'Bad value of datasetId, should be 1seg or 3seg');
end

if ~exist('matlabResultFolder', 'var') || isempty(matlabResultFolder)
    matlabResultFolder = fullfile( 'data', 'Msrc', ['MatlabResults_', datasetId] );
end
if ~ischar(matlabResultFolder) || ~exist(matlabResultFolder, 'dir')
    error('makePlots_MSRC:badMatlabResultsFolder', 'matlabResultsFolder should be a string with a path to the results');
end

if ~exist('resultImageNames', 'var') || isempty(resultImageNames)
    listOfFiles = dir(fullfile(matlabResultFolder, '*.mat'));
    resultImageNames = cell(length(listOfFiles), 1);
    for iImage = 1 : length(listOfFiles)
        [~, resultImageNames{iImage}, ~] = fileparts(listOfFiles(iImage).name);
    end
end
if ~iscellstr(resultImageNames)
    error('makePlots_MSRC:badResultImageNames', 'resultImageNames should be a cell array of strings');
end
resultImageNames = resultImageNames(:);

%% load results
numResults = length(resultImageNames);
results = cell(numResults, 1);

timePlots_smr = cell(numResults, 1);
timePlots_cwd = cell(numResults, 1);
ticPlots_cmr = cell(numResults, 1);
ticPlots_cwd = cell(numResults, 1);
dualPlots_smr = cell(numResults, 1);
dualPlots_cwd = cell(numResults, 1);
primalPlots_smr = cell(numResults, 1);
primalPlots_cwd = cell(numResults, 1);

for iImage = 1 : numResults
    curFileName = fullfile( matlabResultFolder, [resultImageNames{iImage},'.mat'] );
    if ~exist(curFileName, 'file')
        error('makePlots_MSRC:badResultFiles', ['Result file does not exist: ', curFileName ]);
    end
    results{ iImage } = load( curFileName, 'energy_unaryOnly',  'hammingError_unaryOnly', 'energy_unaryPairwise', 'hammingError_unaryPairwise', ...
        'energy_unaryPairwiseHo', 'hammingError_unaryPairwiseHo', 'timePlot_smr', 'dualPlot_smr', 'primalPlot_smr', 'hammingError_unaryPairwiseHo_smr', ...
        'timePlot_cwd', 'dualPlot_cwd', 'primalPlot_cwd', 'hammingError_unaryPairwiseHo_cwd', 'timeAlphaExp', 'timeSmr', 'timeCwd' );
    
    timePlots_smr{ iImage } = results{ iImage }.timePlot_smr;
    timePlots_cwd{ iImage } = results{ iImage }.timePlot_cwd;
    ticPlots_cmr{ iImage } = 1 : length(results{ iImage }.timePlot_smr);
    ticPlots_cwd{ iImage } = 1 : length(results{ iImage }.timePlot_cwd);
    dualPlots_smr{ iImage } = results{ iImage }.dualPlot_smr;
    dualPlots_cwd{ iImage } = results{ iImage }.dualPlot_cwd;
    primalPlots_smr{ iImage } = results{ iImage }.primalPlot_smr;
    primalPlots_cwd{ iImage } = results{ iImage }.primalPlot_cwd;
    
    energy_unaryOnly( iImage ) = results{ iImage }.energy_unaryOnly;
    energy_unaryPairwise( iImage ) = results{ iImage }.energy_unaryPairwise;
    energy_unaryPairwiseHo( iImage ) = results{ iImage }.energy_unaryPairwiseHo;
    bestPrimal_smr( iImage ) = min( results{ iImage }.primalPlot_smr );
    bestDual_smr( iImage ) = max( results{ iImage }.dualPlot_smr );
    bestPrimal_cwd( iImage ) = min( results{ iImage }.primalPlot_cwd );
    bestDual_cwd( iImage ) = max( results{ iImage }.dualPlot_cwd );
    
    time_smr( iImage ) = results{ iImage }.timeSmr;
    time_cwd( iImage ) = results{ iImage }.timeCwd;
    time_alphaExp( iImage ) = results{ iImage }.timeAlphaExp;
end

% subtract the zero level
for iImage = 1 : numResults
    zeroLevel = max( dualPlots_smr{iImage}(~isnan(dualPlots_smr{iImage})));
    
    dualPlots_smr{iImage} = dualPlots_smr{iImage} - zeroLevel;
    primalPlots_smr{iImage} = primalPlots_smr{iImage} - zeroLevel;
    dualPlots_cwd{iImage} = dualPlots_cwd{iImage} - zeroLevel;
    primalPlots_cwd{iImage} = primalPlots_cwd{iImage} - zeroLevel;
    
    energy_unaryOnly( iImage ) = energy_unaryOnly( iImage ) - zeroLevel;
    energy_unaryPairwise( iImage ) = energy_unaryPairwise( iImage ) - zeroLevel;
    energy_unaryPairwiseHo( iImage ) = energy_unaryPairwiseHo( iImage )  - zeroLevel;
    bestPrimal_smr( iImage ) = bestPrimal_smr( iImage ) - zeroLevel;
    bestDual_smr( iImage ) = bestDual_smr( iImage ) - zeroLevel;
    bestPrimal_cwd( iImage ) = bestPrimal_cwd( iImage ) - zeroLevel;
    bestDual_cwd( iImage ) = bestDual_cwd( iImage ) - zeroLevel;
end

% comparison against alpha-expansion
curEps = 1e-6;
curMask = true(numResults, 1);
fprintf('SMR certificate: %d out of %d\n', sum(abs(bestPrimal_smr(curMask) - bestDual_smr(curMask)) < curEps), sum(curMask) );
fprintf('Alpha-expansion got same energy: %d out of %d\n', sum(abs(energy_unaryPairwiseHo(curMask) - bestDual_smr(curMask)) < curEps), sum(curMask) );

curMask = curMask & ( (abs(bestPrimal_smr(:) - bestDual_smr(:)) > curEps) | (abs(energy_unaryPairwiseHo(:) - bestDual_smr(:) ) > curEps) );
[energyMedian_smr, energyLowerQuartile_smr, energyUpperQuartile_smr] = getQuartile(bestPrimal_smr(curMask));
[energyMedian_cwd, energyLowerQuartile_cwd, energyUpperQuartile_cwd] = getQuartile(bestPrimal_cwd(curMask));
[energyMedian_alpha, energyLowerQuartile_alpha, energyUpperQuartile_alpha] = getQuartile(energy_unaryPairwiseHo(curMask));

fprintf('SMR quartiles: %f - lower, %f - median, %f - upper\n', energyLowerQuartile_smr, energyMedian_smr, energyUpperQuartile_smr);
fprintf('Alpha-expansion quartiles: %f - lower, %f - median, %f - upper\n', energyLowerQuartile_alpha, energyMedian_alpha, energyUpperQuartile_alpha);
fprintf('CWD quartiles: %f - lower, %f - median, %f - upper\n', energyLowerQuartile_cwd, energyMedian_cwd, energyUpperQuartile_cwd);

noGapMask = curMask & ( abs(bestPrimal_smr(:) - bestDual_smr(:)) < curEps );
gapMask = curMask & ( abs(bestPrimal_smr(:) - bestDual_smr(:)) >= curEps );

fprintf('SMR median runtime (no gap): %f\n', median( time_smr( noGapMask ) ) );
fprintf('SMR median runtime (gap): %f\n', median( time_smr( gapMask ) ) );
fprintf('CWD median runtime: %f\n', median(time_cwd));
fprintf('alpha-expansion median runtime: %f\n', median(time_alphaExp));

% make the plot
figure;
hold on;

alpha = 0.25;
numTics = 500;
plotSteps = (1 : 1 : numTics);

[xDualSmd, yDualSmd, yDualSmdBoundLower, yDualSmdBoundUpper ] = averagePlotsRobust( ticPlots_cmr, dualPlots_smr, alpha, 0, numTics );
[xDualCwd, yDualCwd, yDualCwdBoundLower, yDualCwdBoundUpper ] = averagePlotsRobust( ticPlots_cwd, dualPlots_cwd, alpha, 0, numTics );
[xPrimalSmr, yPrimalSmr, yPrimalSmrBoundLower, yPrimalSmrBoundUpper ] = averagePlotsRobust( ticPlots_cmr, primalPlots_smr, alpha, 0, numTics );
[xPrimalCwd, yPrimalCwd, yPrimalCwdBoundLower, yPrimalCwdBoundUpper ] = averagePlotsRobust( ticPlots_cwd, primalPlots_cwd, alpha, 0, numTics );

lineWidthMain = 1.5;
lineWidthQuantiles = 1.0;

plot(xDualSmd(plotSteps), yDualSmd(plotSteps), 'r', 'LineWidth', lineWidthMain);
plot(xDualSmd(plotSteps), yDualSmdBoundLower(plotSteps), 'r--', 'LineWidth', lineWidthQuantiles);
plot(xDualSmd(plotSteps), yDualSmdBoundUpper(plotSteps), 'r--', 'LineWidth', lineWidthQuantiles);

plot(xPrimalSmr(plotSteps), yPrimalSmr(plotSteps), 'r', 'LineWidth', lineWidthMain);
plot(xPrimalSmr(plotSteps), yPrimalSmrBoundLower(plotSteps), 'r--', 'LineWidth', lineWidthQuantiles);
plot(xPrimalSmr(plotSteps), yPrimalSmrBoundUpper(plotSteps), 'r--', 'LineWidth', lineWidthQuantiles);

ticSteps = 11 : 20 : numTics;
curMarkerSize = 5;
plot(xPrimalSmr(ticSteps), yPrimalSmr(ticSteps), 'r', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerSize', curMarkerSize);
plot(xPrimalSmr(ticSteps), yPrimalSmrBoundLower(ticSteps), 'r', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerSize', curMarkerSize);
plot(xPrimalSmr(ticSteps), yPrimalSmrBoundUpper(ticSteps), 'r', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerSize', curMarkerSize);

    
plot(xDualCwd(plotSteps), yDualCwd(plotSteps), 'b', 'LineWidth', lineWidthMain);
plot(xDualCwd(plotSteps), yDualCwdBoundLower(plotSteps), 'b--', 'LineWidth', lineWidthQuantiles);
plot(xDualCwd(plotSteps), yDualCwdBoundUpper(plotSteps), 'b--', 'LineWidth', lineWidthQuantiles);

plot(xPrimalCwd(plotSteps), yPrimalCwd(plotSteps), 'b', 'LineWidth', lineWidthMain);
plot(xPrimalCwd(plotSteps), yPrimalCwdBoundLower(plotSteps), 'b--', 'LineWidth', lineWidthQuantiles);
plot(xPrimalCwd(plotSteps), yPrimalCwdBoundUpper(plotSteps), 'b--', 'LineWidth', lineWidthQuantiles);

ticSteps = 11 : 20 : numTics;
curMarkerSize = 5;
plot(xPrimalCwd(ticSteps), yPrimalCwd(ticSteps), 'b', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerSize', curMarkerSize);
plot(xPrimalCwd(ticSteps), yPrimalCwdBoundLower(ticSteps), 'b', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerSize', curMarkerSize);
plot(xPrimalCwd(ticSteps), yPrimalCwdBoundUpper(ticSteps), 'b', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerSize', curMarkerSize);

set(gca, 'xlim', [0,300]);
set(gca, 'ylim', [-900, 10000]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 10);

title(['Comparison of SMR vs CWD, dataset: ', datasetId])
ylabel('energies and lower bounds');
xlabel('number of oracle calls');

hold off;

end
