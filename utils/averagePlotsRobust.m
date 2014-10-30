function [ xMean, yMean, yBoundLower, yBoundUpper ] = averagePlotsRobust( xArgs, yArgs, alpha, monotoneFlag, numTics )
%averagePlotsRobust averages plots (with resamling) in a robust way.
% xArgs, yArgs - cell arrays of x- and y-coordinates of the plots
% alpha% of top and alpha% of worst points are not taken into account
% monotoneFlag: 0 do nothing, 1 - make increasing, -1 - make decreasing (default: 0)
% numTics: number of tics to resample the plots
% 
% 	Anton Osokin (firstname.lastname@gmail.com),  14.10.2014

if ~exist( 'monotoneFlag', 'var')
    monotoneFlag = 0;
end
if ~exist( 'numTics', 'var' )
    numTics = 100;
end

xArgs = xArgs(:);
yArgs = yArgs(:);

numPlots = length( xArgs );

% convert to cell array
if ~iscell(xArgs)
    temp = xArgs;
    xArgs = cell( length( yArgs ), 1 );
    for iPlot = 1 : numPlots
        xArgs{iPlot} = temp;
    end
end    

% find min and max args
argmin = inf;
argmax = -inf;

for iPlot = 1 : numPlots
    argumentNanMask = isnan( xArgs{iPlot} );
    
    xArgs{iPlot}( argumentNanMask ) = [];
    yArgs{iPlot}( argumentNanMask ) = [];
    
    if monotoneFlag == 1
        for iTic = 2 : length(yArgs{iPlot})
            yArgs{iPlot}(iTic) = max(yArgs{iPlot}(iTic), yArgs{iPlot}(iTic - 1));
        end
    end

    if monotoneFlag == -1
        for iTic = 2 : length(yArgs{iPlot})
            yArgs{iPlot}(iTic) = min(yArgs{iPlot}(iTic), yArgs{iPlot}(iTic - 1));
        end
    end

    argmin = min(argmin, min( xArgs{ iPlot } ));
    argmax = max(argmax, max( xArgs{ iPlot } ));
end



% new arguments
xMean = argmin : (argmax - argmin) / (numTics - 1) : argmax;

yArgsNew = nan( numPlots, numTics );

for iPlot = 1 : numPlots
    lastNonNanId = find(~isnan(yArgs{iPlot}), 1, 'last');
    firstNonNanId = find(~isnan(yArgs{iPlot}), 1, 'first');
    yArgs{iPlot}( isnan(yArgs{iPlot}(:)) & lastNonNanId < (1 : length(yArgs{iPlot}))' ) = yArgs{iPlot}( lastNonNanId );
    yArgs{iPlot}( isnan(yArgs{iPlot}(:)) & firstNonNanId >(1 : length(yArgs{iPlot}))' ) = yArgs{iPlot}( firstNonNanId );
    
    if length(xArgs{iPlot}) == 1
        xArgs{iPlot} = [xArgs{iPlot}; xArgs{iPlot} + 1];
        yArgs{iPlot} = [yArgs{iPlot}; yArgs{iPlot}];
    end
    temp = interp1( xArgs{iPlot}, yArgs{iPlot}, xMean, 'nearest', nan );
    
    nonNanId = find( ~isnan( temp ), 1, 'first' );
    nanId = find( isnan( temp ) );
    
    
    temp( nanId( nanId < nonNanId ) ) = yArgs{iPlot}( firstNonNanId );
    temp( nanId( nanId > nonNanId ) ) = yArgs{iPlot}( lastNonNanId );
    
    yArgsNew(iPlot, :) = temp;
end

yMean = nan(size(xMean));
yBoundLower = nan(size(xMean));
yBoundUpper = nan(size(xMean));

for iTic = 1 : numTics
    curData = yArgsNew(:, iTic);
    
    curData = sort(curData);
    dataSize = length(curData);
    
    cutOffNumber = round(dataSize * alpha);
    curData(1 : cutOffNumber) = [];
    curData(end - cutOffNumber + 1 : end) = [];
    
    yMean(iTic) = median(curData);
    yBoundUpper(iTic) = curData(end);
    yBoundLower(iTic) = curData(1);
end

end

