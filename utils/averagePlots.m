function [ xMean, yMean, yStd ] = averagePlots( xArgs, yArgs, monotoneFlag )
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

yMean = mean( yArgsNew, 1 ); 
yStd =  std(  yArgsNew, 1 ); 

end

