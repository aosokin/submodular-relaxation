function [ xMean, yMean, yStd ] = averagePlots( xArgs, yArgs, monotoneFlag, numTics, timeInterval )
%averagePlots averages several plots together with resampling
% xArgs, yArgs - cell arrays of x- and y-coordinates of the plots
% monotoneFlag: 0 do nothing, 1 - make increasing, -1 - make decreasing (default: 0)
% numTics - number of tics in the plots (default: 100)
% timeInterval - vector of length 2 with starting and ending points of the X-axis (default: unused)
% 
% 	Anton Osokin (firstname.lastname@gmail.com),  14.10.2014

%% parse the input 
if ~exist( 'monotoneFlag', 'var') || isempty(monotoneFlag)
    monotoneFlag = 0;
end
if ~isnumeric(monotoneFlag) || ~isscalar(monotoneFlag) || ~ismember(monotoneFlag, [-1, 0, 1] )
    error('averagePlots:badMonotoneFlag', 'monotoneFlag should a single number form set [-1, 0, 1]');
end

if ~exist( 'numTics', 'var') || isempty(numTics)
    numTics = 100;
end
if ~isnumeric(numTics) || ~isscalar(numTics) || numTics < 2
    error('averagePlots:badNumTics', 'numTics (the number of tics) should a single positive number, numTics >= 2');
end
numTics = round(numTics);

if ~exist( 'timeInterval', 'var')
    timeInterval = [];
end
if ~isempty(timeInterval) && (~isnumeric(timeInterval) || ~isvector(timeInterval) || length(timeInterval) ~= 2 || any(isnan(timeInterval)) || any(isinf(timeInterval)) || any(timeInterval < 0) || timeInterval(1) >= timeInterval(2))
    error('averagePlots:badTimeInterval', 'timeInterval should be an increasing vector of two non-negative elemets');
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

%% start averaging

% find min and max args
argmin = -inf;
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
    
    argmin = max(argmin, min( xArgs{ iPlot } ));
    argmax = max(argmax, max( xArgs{ iPlot } ));
end

if isempty(timeInterval)
    timeInterval = [argmin, argmax];
end

% new arguments
xMean = timeInterval(1) : (timeInterval(2) - timeInterval(1)) / (numTics - 1) : timeInterval(2);

yArgsNew = nan( numPlots, numTics );

for iPlot = 1 : numPlots
    lastNonNanId = find(~isnan(yArgs{iPlot}), 1, 'last');
    firstNonNanId = find(~isnan(yArgs{iPlot}), 1, 'first');
    yArgs{iPlot}( isnan(yArgs{iPlot}(:)) & lastNonNanId < (1 : length(yArgs{iPlot}))' ) = yArgs{iPlot}( lastNonNanId );
    yArgs{iPlot}( isnan(yArgs{iPlot}(:)) & firstNonNanId >(1 : length(yArgs{iPlot}))' ) = yArgs{iPlot}( firstNonNanId );
    
    temp = interp1( xArgs{iPlot}, yArgs{iPlot}, xMean, 'nearest', nan );
    
    nonNanId = find( ~isnan( temp ), 1, 'first' );
    nanId = find( isnan( temp ) );
    
    
    temp( nanId( nanId < nonNanId ) ) = yArgs{iPlot}( firstNonNanId );
    temp( nanId( nanId > nonNanId ) ) = yArgs{iPlot}( lastNonNanId );
    
    yArgsNew(iPlot, :) = temp;
end

yMean = mean( yArgsNew, 1 );
yStd =  std(  yArgsNew, 1 );

end

