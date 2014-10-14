function makePlots_SMR_stereo(resultFileName)
% makePlots_SMR_stereo plots fig.4b from the journal paper
%   resultFileName - file with the experiment results (default: results_SMR.mat)
%
% Anton Osokin (firstname.lastname@gmail.com),  14.10.2014

setup_SMR;

if ~exist('resultFileName', 'var')
    resultFileName = 'results_SMR.mat';
end
objectIds = [1, 2, 3, 4];

%% load results
load(resultFileName, 'timePlot_trws' ,'dualPlot_trws', 'primalPlot_trws',  ...
               'timePlot_SMR_Hanso', 'dualPlot_SMR_Hanso', 'primalPlot_SMR_Hanso', ...
               'timePlot_SMR_LMBM', 'dualPlot_SMR_LMBM', 'primalPlot_SMR_LMBM', ...
               'timePlot_SMR_bundle', 'dualPlot_SMR_bundle', 'primalPlot_SMR_bundle', ...
               'timePlot_SMR_bundleAggr', 'dualPlot_SMR_bundleAggr', 'primalPlot_SMR_bundleAggr', ...
               'timePlot_SMR_subgradient', 'dualPlot_SMR_subgradient', 'primalPlot_SMR_subgradient', ...
               'timePlot_DDTRW_Hanso', 'dualPlot_DDTRW_Hanso', 'primalPlot_DDTRW_Hanso', ...
               'timePlot_DDTRW_LMBM', 'dualPlot_DDTRW_LMBM', 'primalPlot_DDTRW_LMBM', ...
               'timePlot_DDTRW_bundle', 'dualPlot_DDTRW_bundle', 'primalPlot_DDTRW_bundle', ...
               'timePlot_DDTRW_bundleAggr', 'dualPlot_DDTRW_bundleAggr', 'primalPlot_DDTRW_bundleAggr', ...
               'timePlot_DDTRW_subgradient', 'dualPlot_DDTRW_subgradient', 'primalPlot_DDTRW_subgradient');


%% average plots
numTics = 100;
maxTime = 200;
timeInterval = [0.1, maxTime];

[ticsDualPlotFinal_trws, dualPlotFinal_trws] = averagePlots( timePlot_trws(objectIds), dualPlot_trws(objectIds), 1, numTics, timeInterval );
[ticsPrimalPlotFinal_trws, primalPlotFinal_trws] = averagePlots( timePlot_trws(objectIds), primalPlot_trws(objectIds), -1, numTics, timeInterval );

[ticsDualPlotFinal_SMR_subgradient, dualPlotFinal_SMR_subgradient] = averagePlots( timePlot_SMR_subgradient(objectIds), dualPlot_SMR_subgradient(objectIds), 1, numTics, timeInterval );
[ticsPrimalPlotFinal_SMR_subgradient, primalPlotFinal_SMR_subgradient] = averagePlots( timePlot_SMR_subgradient(objectIds), primalPlot_SMR_subgradient(objectIds), -1, numTics, timeInterval );

[ticsDualPlotFinal_SMR_bundleAggr, dualPlotFinal_SMR_bundleAggr] = averagePlots( timePlot_SMR_bundleAggr(objectIds), dualPlot_SMR_bundleAggr(objectIds), 1, numTics, timeInterval );
[ticsPrimalPlotFinal_SMR_bundleAggr, primalPlotFinal_SMR_bundleAggr] = averagePlots( timePlot_SMR_bundleAggr(objectIds), primalPlot_SMR_bundleAggr(objectIds), -1, numTics, timeInterval );

[ticsDualPlotFinal_SMR_bundle, dualPlotFinal_SMR_bundle] = averagePlots( timePlot_SMR_bundle(objectIds), dualPlot_SMR_bundle(objectIds), 1, numTics, timeInterval );
[ticsPrimalPlotFinal_SMR_bundle, primalPlotFinal_SMR_bundle] = averagePlots( timePlot_SMR_bundle(objectIds), primalPlot_SMR_bundle(objectIds), -1, numTics, timeInterval );

[ticsDualPlotFinal_SMR_LMBM, dualPlotFinal_SMR_LMBM] = averagePlots( timePlot_SMR_LMBM(objectIds), dualPlot_SMR_LMBM(objectIds), 1, numTics, timeInterval );
[ticsPrimalPlotFinal_SMR_LMBM, primalPlotFinal_SMR_LMBM] = averagePlots( timePlot_SMR_LMBM(objectIds), primalPlot_SMR_LMBM(objectIds), -1, numTics, timeInterval );

[ticsDualPlotFinal_SMR_Hanso, dualPlotFinal_SMR_Hanso] = averagePlots( timePlot_SMR_Hanso(objectIds), dualPlot_SMR_Hanso(objectIds), 1, numTics, timeInterval );
[ticsPrimalPlotFinal_SMR_Hanso, primalPlotFinal_SMR_Hanso] = averagePlots( timePlot_SMR_Hanso(objectIds), primalPlot_SMR_Hanso(objectIds), -1, numTics, timeInterval );

[ticsDualPlotFinal_DDTRW_subgradient, dualPlotFinal_DDTRW_subgradient] = averagePlots( timePlot_DDTRW_subgradient(objectIds), dualPlot_DDTRW_subgradient(objectIds), 1, numTics, timeInterval );
[ticsPrimalPlotFinal_DDTRW_subgradient, primalPlotFinal_DDTRW_subgradient] = averagePlots( timePlot_DDTRW_subgradient(objectIds), primalPlot_DDTRW_subgradient(objectIds), -1, numTics, timeInterval );

[ticsDualPlotFinal_DDTRW_bundleAggr, dualPlotFinal_DDTRW_bundleAggr] = averagePlots( timePlot_DDTRW_bundleAggr(objectIds), dualPlot_DDTRW_bundleAggr(objectIds), 1, numTics, timeInterval );
[ticsPrimalPlotFinal_DDTRW_bundleAggr, primalPlotFinal_DDTRW_bundleAggr] = averagePlots( timePlot_DDTRW_bundleAggr(objectIds), primalPlot_DDTRW_bundleAggr(objectIds), -1, numTics, timeInterval );

[ticsDualPlotFinal_DDTRW_bundle, dualPlotFinal_DDTRW_bundle] = averagePlots( timePlot_DDTRW_bundle(objectIds), dualPlot_DDTRW_bundle(objectIds), 1, numTics, timeInterval );
[ticsPrimalPlotFinal_DDTRW_bundle, primalPlotFinal_DDTRW_bundle] = averagePlots( timePlot_DDTRW_bundle(objectIds), primalPlot_DDTRW_bundle(objectIds), -1, numTics, timeInterval );

[ticsDualPlotFinal_DDTRW_LMBM, dualPlotFinal_DDTRW_LMBM] = averagePlots( timePlot_DDTRW_LMBM(objectIds), dualPlot_DDTRW_LMBM(objectIds), 1, numTics, timeInterval );
[ticsPrimalPlotFinal_DDTRW_LMBM, primalPlotFinal_DDTRW_LMBM] = averagePlots( timePlot_DDTRW_LMBM(objectIds), primalPlot_DDTRW_LMBM(objectIds), -1, numTics, timeInterval );

[ticsDualPlotFinal_DDTRW_Hanso, dualPlotFinal_DDTRW_Hanso] = averagePlots( timePlot_DDTRW_Hanso(objectIds), dualPlot_DDTRW_Hanso(objectIds), 1, numTics, timeInterval );
[ticsPrimalPlotFinal_DDTRW_Hanso, primalPlotFinal_DDTRW_Hanso] = averagePlots( timePlot_DDTRW_Hanso(objectIds), primalPlot_DDTRW_Hanso(objectIds), -1, numTics, timeInterval );


%% Dual plots
figure(3);
clf;
hold on


lineWidth = 1.5;
curMarkerSize = 7;

legendNames = cell(0, 0);
markerPlotX = cell(0, 0);
markerPlotY = cell(0, 0);
markerType = cell(0,0);
markerColor = cell(0, 0);
dataPlotX = cell(0, 0);
dataPlotY = cell(0, 0);
dataPlotColor = cell(0, 0);
curPlotIndex = 0;

% TRW-S
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsDualPlotFinal_trws;
dataPlotY{curPlotIndex} = dualPlotFinal_trws;
legendNames{curPlotIndex} = 'TRW-S';
curSteps = findClosest(dataPlotX{curPlotIndex}, [2, 5 : 60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex} = 'p';
markerColor{curPlotIndex} = 'k';
dataPlotColor{curPlotIndex} = 'k';

% SMR bundle
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsDualPlotFinal_SMR_bundle;
dataPlotY{curPlotIndex} = dualPlotFinal_SMR_bundle;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'SMR bundle';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 80, 10 : 60: maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'o';
markerColor{curPlotIndex} = 'r';
dataPlotColor{curPlotIndex} = 'r';

% SMR aggr. bundle
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsDualPlotFinal_SMR_bundleAggr;
dataPlotY{curPlotIndex} = dualPlotFinal_SMR_bundleAggr;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'SMR aggr. bundle';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 20 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 's';
markerColor{curPlotIndex} = 'r';
dataPlotColor{curPlotIndex} = 'r';

% SMR subgradient
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsDualPlotFinal_SMR_subgradient;
dataPlotY{curPlotIndex} = dualPlotFinal_SMR_subgradient;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'SMR subgradient';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 40,30 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= '^';
markerColor{curPlotIndex} = 'r';
dataPlotColor{curPlotIndex} = 'r';

% SMR Hanso
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsDualPlotFinal_SMR_Hanso;
dataPlotY{curPlotIndex} = dualPlotFinal_SMR_Hanso;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'SMR L-BFGS';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 45, 40 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'v';
markerColor{curPlotIndex} = 'r';
dataPlotColor{curPlotIndex} = 'r';

% SMR LMBM
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsDualPlotFinal_SMR_LMBM;
dataPlotY{curPlotIndex} = dualPlotFinal_SMR_LMBM;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'SMR LMBM';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 50 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'd';
markerColor{curPlotIndex} = 'r';
dataPlotColor{curPlotIndex} = 'r';

% DD TRW bundle
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsDualPlotFinal_DDTRW_bundle;
dataPlotY{curPlotIndex} =  dualPlotFinal_DDTRW_bundle;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'DD TRW bundle';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 15 : 60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'o';
markerColor{curPlotIndex} = 'b';
dataPlotColor{curPlotIndex} = 'b';

% DD TRW aggr. bundle
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsDualPlotFinal_DDTRW_bundleAggr;
dataPlotY{curPlotIndex} = dualPlotFinal_DDTRW_bundleAggr;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'DD TRW aggr. bundle';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 40, 25 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex} = 's';
markerColor{curPlotIndex} = 'b';
dataPlotColor{curPlotIndex} = 'b';

% DD TRW subgradient
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsDualPlotFinal_DDTRW_subgradient;
dataPlotY{curPlotIndex} = dualPlotFinal_DDTRW_subgradient;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'DD TRW subgradient';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 35 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= '^';
markerColor{curPlotIndex} = 'b';
dataPlotColor{curPlotIndex} = 'b';

% DD TRW Hanso
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsDualPlotFinal_DDTRW_Hanso;
dataPlotY{curPlotIndex} =  dualPlotFinal_DDTRW_Hanso;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'DD TRW L-BFGS';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 45 :  50 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'v';
markerColor{curPlotIndex} = 'b';
dataPlotColor{curPlotIndex} = 'b';

% DD TRW LMBM
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsDualPlotFinal_DDTRW_LMBM;
dataPlotY{curPlotIndex} = dualPlotFinal_DDTRW_LMBM;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'DD TRW LMBM';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 55 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'd';
markerColor{curPlotIndex} = 'b';
dataPlotColor{curPlotIndex} = 'b';


hold on;
% plot lines first 
for iPlot = 1 : 1 : curPlotIndex
    plot(dataPlotX{iPlot}, dataPlotY{iPlot}, 'Color', dataPlotColor{iPlot},  'LineStyle', '-', 'LineWidth', lineWidth);    
end

% plot the markers
for iPlot = 1 : 1 : curPlotIndex
    plot(markerPlotX{iPlot}, markerPlotY{iPlot}, 'Color', dataPlotColor{iPlot},  'LineStyle', 'none', 'Marker', markerType{iPlot}, 'MarkerFaceColor', markerColor{iPlot}, 'MarkerSize', curMarkerSize, 'MarkerEdgeColor', 'k');    
end

% generate dummy plots for legends
legendHandles = zeros(curPlotIndex,1);
for iPlot = 1 : 1 : curPlotIndex
    legendHandles(iPlot) = plot(0, 0, 'Marker', markerType{iPlot}, 'MarkerFaceColor', markerColor{iPlot}, 'MarkerSize', curMarkerSize, 'MarkerEdgeColor', 'k', 'LineStyle', '-', 'LineWidth', lineWidth, 'Color', dataPlotColor{iPlot});    
end
hold off;

% adjust the plot
legend( legendHandles, legendNames );
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 10);
xlim([0 maxTime]);
ylim([97, 100.2]);
ylabel('lower bounds');
xlabel('time (seconds)');
title('SMR and DD-TRW lower bounds for stereo instances');


%% Primal plots
figure(4);
clf;
hold on


lineWidth = 1.5;
curMarkerSize = 7;

legendNames = cell(0, 0);
markerPlotX = cell(0, 0);
markerPlotY = cell(0, 0);
markerType = cell(0,0);
markerColor = cell(0, 0);
dataPlotX = cell(0, 0);
dataPlotY = cell(0, 0);
dataPlotColor = cell(0, 0);
curPlotIndex = 0;

% TRW-S
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsPrimalPlotFinal_trws;
dataPlotY{curPlotIndex} = primalPlotFinal_trws;
legendNames{curPlotIndex} = 'TRW-S';
curSteps = findClosest(dataPlotX{curPlotIndex}, [2, 5 : 60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex} = 'p';
markerColor{curPlotIndex} = 'k';
dataPlotColor{curPlotIndex} = 'k';

% SMR bundle
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsPrimalPlotFinal_SMR_bundle;
dataPlotY{curPlotIndex} = primalPlotFinal_SMR_bundle;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'SMR bundle';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 80, 10 : 60: maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'o';
markerColor{curPlotIndex} = 'r';
dataPlotColor{curPlotIndex} = 'r';

% SMR aggr. bundle
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsPrimalPlotFinal_SMR_bundleAggr;
dataPlotY{curPlotIndex} = primalPlotFinal_SMR_bundleAggr;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'SMR aggr. bundle';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 20 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 's';
markerColor{curPlotIndex} = 'r';
dataPlotColor{curPlotIndex} = 'r';

% SMR subgradient
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsPrimalPlotFinal_SMR_subgradient;
dataPlotY{curPlotIndex} = primalPlotFinal_SMR_subgradient;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'SMR subgradient';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 40,30 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= '^';
markerColor{curPlotIndex} = 'r';
dataPlotColor{curPlotIndex} = 'r';

% SMR Hanso
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsPrimalPlotFinal_SMR_Hanso;
dataPlotY{curPlotIndex} = primalPlotFinal_SMR_Hanso;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'SMR L-BFGS';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 45, 40 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'v';
markerColor{curPlotIndex} = 'r';
dataPlotColor{curPlotIndex} = 'r';

% SMR LMBM
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsPrimalPlotFinal_SMR_LMBM;
dataPlotY{curPlotIndex} = primalPlotFinal_SMR_LMBM;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'SMR LMBM';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 50 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'd';
markerColor{curPlotIndex} = 'r';
dataPlotColor{curPlotIndex} = 'r';

% DD TRW bundle
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsPrimalPlotFinal_DDTRW_bundle;
dataPlotY{curPlotIndex} =  primalPlotFinal_DDTRW_bundle;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'DD TRW bundle';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 15 : 60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'o';
markerColor{curPlotIndex} = 'b';
dataPlotColor{curPlotIndex} = 'b';

% DD TRW aggr. bundle
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsPrimalPlotFinal_DDTRW_bundleAggr;
dataPlotY{curPlotIndex} = primalPlotFinal_DDTRW_bundleAggr;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'DD TRW aggr. bundle';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 40, 25 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex} = 's';
markerColor{curPlotIndex} = 'b';
dataPlotColor{curPlotIndex} = 'b';

% DD TRW subgradient
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsPrimalPlotFinal_DDTRW_subgradient;
dataPlotY{curPlotIndex} = primalPlotFinal_DDTRW_subgradient;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'DD TRW subgradient';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 35 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= '^';
markerColor{curPlotIndex} = 'b';
dataPlotColor{curPlotIndex} = 'b';

% DD TRW Hanso
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsPrimalPlotFinal_DDTRW_Hanso;
dataPlotY{curPlotIndex} =  primalPlotFinal_DDTRW_Hanso;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'DD TRW L-BFGS';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 45 :  50 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'v';
markerColor{curPlotIndex} = 'b';
dataPlotColor{curPlotIndex} = 'b';

% DD TRW LMBM
curPlotIndex = curPlotIndex + 1;
dataPlotX{curPlotIndex} = ticsPrimalPlotFinal_DDTRW_LMBM;
dataPlotY{curPlotIndex} = primalPlotFinal_DDTRW_LMBM;
if max(dataPlotX{curPlotIndex}) < maxTime
    dataPlotX{curPlotIndex} = [dataPlotX{curPlotIndex}, max(dataPlotX{curPlotIndex}) : 0.1: maxTime];
    dataPlotY{curPlotIndex} = [dataPlotY{curPlotIndex}, max(dataPlotY{curPlotIndex}) * ones(1, size(dataPlotX{curPlotIndex}, 2) - size(dataPlotY{curPlotIndex},2))];
end
legendNames{curPlotIndex} = 'DD TRW LMBM';
curSteps = findClosest(dataPlotX{curPlotIndex}, [.8, 55 :  60 : maxTime]);
markerPlotX{curPlotIndex} = dataPlotX{curPlotIndex}(curSteps);
markerPlotY{curPlotIndex} = dataPlotY{curPlotIndex}(curSteps);
markerType{curPlotIndex}= 'd';
markerColor{curPlotIndex} = 'b';
dataPlotColor{curPlotIndex} = 'b';


hold on;
% plot lines first 
for iPlot = 1 : 1 : curPlotIndex
    plot(dataPlotX{iPlot}, dataPlotY{iPlot}, 'Color', dataPlotColor{iPlot},  'LineStyle', '-', 'LineWidth', lineWidth);    
end

% plot the markers
for iPlot = 1 : 1 : curPlotIndex
    plot(markerPlotX{iPlot}, markerPlotY{iPlot}, 'Color', dataPlotColor{iPlot},  'LineStyle', 'none', 'Marker', markerType{iPlot}, 'MarkerFaceColor', markerColor{iPlot}, 'MarkerSize', curMarkerSize, 'MarkerEdgeColor', 'k');    
end

% generate dummy plots for legends
legendHandles = zeros(curPlotIndex,1);
for iPlot = 1 : 1 : curPlotIndex
    legendHandles(iPlot) = plot(0, 0, 'Marker', markerType{iPlot}, 'MarkerFaceColor', markerColor{iPlot}, 'MarkerSize', curMarkerSize, 'MarkerEdgeColor', 'k', 'LineStyle', '-', 'LineWidth', lineWidth, 'Color', dataPlotColor{iPlot});    
end
hold off;

% adjust the plot
legend( legendHandles, legendNames );
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 10);
xlim([0 maxTime]);
ylim([99.8, 150]);
ylabel('energies');
xlabel('time (seconds)');
title('SMR and DD-TRW energies for stereo instances');
