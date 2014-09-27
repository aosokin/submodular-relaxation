function [ xMean, yMean, yBoundLower, yBoundUpper ] = averagePlotsRobust( xArgs, yArgs, alpha, monotoneFlag )
% alpha% of top and alpha% of worst points are not taken into account
% monotoneFlag: 0 do nothing, 1 - make increasing, -1 - make decreasing

if ~exist( 'monotoneFlag', 'var')
    monotoneFlag = 0;
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
numTics = 500;
xMean = argmin : (argmax - argmin) / (numTics - 1) : argmax;

yArgsNew = nan( numPlots, numTics );

for iPlot = 1 : numPlots
    lastNonNan = yArgs{iPlot}( find(~isnan(yArgs{iPlot}), 1, 'last') );
    yArgs{iPlot}( isnan(yArgs{iPlot}) ) = lastNonNan;
    
    temp = interp1( xArgs{iPlot}, yArgs{iPlot}, xMean, 'nearest', nan );
    
    nonNanId = find( ~isnan( temp ), 1, 'first' );
    nanId = find( isnan( temp ) );
    
    
    temp( nanId( nanId < nonNanId ) ) = temp( find( ~isnan( temp ), 1, 'first' ) );
    temp( nanId( nanId > nonNanId ) ) = temp( find( ~isnan( temp ), 1, 'last' ) );
    
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

% yMean = median( yArgsNew, 1 ); 
% yStd =  median(abs(bsxfun(@minus, yArgsNew, yMean)), 1);

end

