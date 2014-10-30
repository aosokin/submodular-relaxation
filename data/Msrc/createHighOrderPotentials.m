function [ finalHoId, finalHoP ] = createHighOrderPotentials( curImage, segments, options )
% readMsrcGroundTruth is a helper function to work with the MSRC-21 dataset 
%
% See createDataset_MSRC for instructions
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

%% check data
curImage = im2double(curImage);

if ~iscell(segments)
    segments = {segments};
end
numSegmentations = length( segments );

imHeight = size( curImage, 1 );
imWidth = size( curImage, 2 );
for iSegmentation = 1 : numSegmentations
    if ~isnumeric( segments{iSegmentation} ) || any( size(segments{ iSegmentation }) ~= [imHeight, imWidth])
        error('createHighOrderPotentials:badSegmentation', ['Segmentation ', num2str(iSegmentation), ' is not compatible with the image']);
    end
    segments{ iSegmentation } = double(segments{ iSegmentation }) - min(segments{ iSegmentation }(:),[],1) + 1;
end

%% parameters
theta_alpha = options.theta_alpha; % 0.8;
theta_hp = options.theta_hp; % 0.2;
theta_hv = options.theta_hv; % 0.5;
theta_hb = options.theta_hb; % 12;
Q_ratio = options.Q_ratio; % 0.1;

%% create potentialss
finalHoId = cell(0, 1);
finalHoP = zeros(0, 2);
for iSegmentation = 1 : numSegmentations
    % superpixels
    curSeg = segments{iSegmentation};
    numSegm = max( curSeg(:) );
    hoId = cell( numSegm, 1 );
    hoP = nan( numSegm, 2 );
    for iSuperpixel = 1 : numSegm
        curMask = (curSeg == iSuperpixel);
        hoId{iSuperpixel} = find( curMask );
        
        numPixels = length( hoId{iSuperpixel} );
        
        curPoints = nan( numPixels, size(curImage, 3));
        for iColor = 1 : size(curImage, 3)
            tmp = curImage(:, :, iColor);
            curPoints(:, iColor) = double( tmp( curMask ) );
        end
        
        mu = mean(curPoints);
        
        coeff = norm( sum( bsxfun( @minus, curPoints, mu ) .^ 2, 1) ) / numPixels;
        
        curG = exp( -theta_hb * coeff );
        
        gammaMax = numPixels ^ theta_alpha * ( theta_hp + theta_hv * curG );
        
        Q = Q_ratio * numPixels;
        
        hoP( iSuperpixel, 1 ) = gammaMax;
        hoP( iSuperpixel, 2 ) = Q;
    end
    
    finalHoId = [ finalHoId; hoId ];
    finalHoP = [ finalHoP; hoP ];
end

end

