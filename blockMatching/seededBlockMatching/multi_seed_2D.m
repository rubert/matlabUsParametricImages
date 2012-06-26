function [dpY, dpX, quality, geom] = multi_seed_2D(pre, post, params)
%%
%I.  DESCRIPTION
%This code is an implementation of the seeds algorithm proposed by Lujie
%Chen et al in A quality guided displacement tracking algorithm for
%ultrasonic elasticity imaging.  This algorithm performs 2-D
%motion tracking axially and laterally.  The quality metric is the normalized
%cross-correlation, and the sub-sample interpolation of the cross
%correlation is done through a parabolic fit, in the axial and lateral
%directions separately.
%
%
%
%II.  OUTPUT
%dpY - the displacements in the y direction, units are pixels
%dpX - displacements in x direction, units are pixels (A-lines)
%quality - The value of the normalized CC associated with each Dp estimate
%geom - a structure containing the following fields:
%        stepY -  spacing between displacement estimates (pixels)
%        startY - location of first displacemetn estimate relative to RF data
%        stopY - location of last displacemetn estimate relative to RF data
%        startX - location of first displacemetn estimate relative to RF data
%        stopX - location of last displacemetn estimate relative to RF data
%
%
%
%
%III.  INPUT
%pre - matrix of RF data before compression
%post - matrix of RF data after compression
%params - a structure I created to keep track of motion tracking
%         paramaters, it must contain the following fields:
%
%         rangeY :  The search range in the axial direction in pixels, used
%         when calculating CC for seeds
%
%         smallRangeY :  The smaller search range in the axial direction
%
%         rangeX :  The search range in the lateral direction in pixels, used
%         when calculating CC for seeds
%
%         smallRangeX :  The smaller search range in the axial direction
%
%         overlap :  The overlap between adjacent blocks axially, a
%         percentage expressed as a number between 0 and 1
%
%         windowY :  The block matching kernel size in the axial
%         direction, THIS NUMBER MUST BE AN ODD INTEGER
%
%         windowX :  The block matching kernel size in the lateral
%         direction, THIS NUMBER MUST BE AN ODD INTEGER
%
%
%
%
%IV.DEPENDENCIES
%To perform the normalized cross correlation rapidly, this script
%relies on a mex file for normalized cross correlation obtainable from:
%http://www.cs.ubc.ca/~deaton/remarks_ncc.html
%The source code for the NCC is taken from the openCV library




%ADJUST PARAMETERS
rangeY = params.rangeY;
halfY = (params.windowY - 1)/2 ;
halfX = (params.windowX - 1)/2;
rangeX = params.rangeX;

smallRangeY = params.smallRangeY;
smallRangeX = params.smallRangeX;


[rfPoints, numLines] = size(pre);




%work out the number of displacement points that will fit
startY = 1 + rangeY + halfY + smallRangeY;
stepY = round(params.windowY*(1-params.overlap));
upperBound = rfPoints - halfY - 1 - rangeY - smallRangeY;
numY = floor( (upperBound - startY) /stepY );

startX = 1 + rangeX + halfX + smallRangeX;
stopX = numLines - 1 - rangeX - halfX - smallRangeX;
numX = 1 + stopX - startX;

%calculate window centers
kernCenterY = zeros(numY, 1);
kernCenterY(1) = startY;

for yy = 2:numY
    
    kernCenterY(yy) = kernCenterY(yy-1) + stepY;
    
end

stopY = kernCenterY(end);

kernCenterX = zeros(numX, 1);
kernCenterX(1) = startX;


for xx = 2:numX
    kernCenterX(xx) = kernCenterX(xx-1) + 1;
end


totPoints = numY*numX;

%the seed step is to chunk the image up into squares
%having a certain percent of the image area
ySeedLength = round( sqrt(.1)*numY );
xSeedLength = round( sqrt(.1)*numX );

thresh = xSeedLength*ySeedLength;

%making sure all the initial seeds are interior points
%to simplify code
seedsY = 2:ySeedLength:numY-1;
seedsX = 2:xSeedLength:numX-1;
totSeeds = length(seedsY)*length(seedsX);


%region set will keep track of which seed a group of points belongs
%to

regionSet = zeros(numY,numX);

dpY = zeros(numY, numX);
dpX = zeros(numY, numX);
quality = -7*ones(numY, numX);
tempQuality = -7*ones(numY, numX);
processed = false(totPoints, 1);


iniDpY = zeros(totPoints,1);
iniDpX = zeros(totPoints,1);

region = 1;

seedPoints = zeros(length(seedsY)*length(seedsX),1 );




%PERFORM CROSS CORRELATION FOR SEEDS
ind = 1;

for yy = 1: length(seedsY)
    for xx = 1:length(seedsX)
        
        
        yInd = seedsY(yy);
        xInd = seedsX(xx);
        pt = sub2ind( [numY, numX], yInd, xInd);
        seedPoints(ind) = pt;
        
        CC = 5 + normxcorr2_mex...
            (pre( kernCenterY(yInd) - halfY : kernCenterY(yInd) + halfY, kernCenterX(xInd) - halfX: kernCenterX(xInd) + halfX )...
            , post(kernCenterY(yInd) - halfY - rangeY : kernCenterY(yInd)  + halfY + rangeY , kernCenterX(xInd) - halfX - rangeX : kernCenterX(xInd) + halfX + rangeX ), 'valid' );
        
        
        [quality(pt), tempInd] = max(CC(:));
        tempQuality(pt) = -7;
        
        
        [dpYInd, dpXInd] = ind2sub(size(CC), tempInd);
        
        if dpYInd > 1 && dpYInd < 2*rangeY + 1
            deltaY = subSampleFit(CC(dpYInd - 1:dpYInd + 1, dpXInd) );
        else
            deltaY = 0;
        end
        
        if dpXInd > 1 && dpXInd < 2*rangeX + 1
            deltaX = subSampleFit(CC(dpYInd, dpXInd - 1: dpXInd + 1) );
        else
            deltaX = 0;
        end
        
        
        
        dpY(pt) = dpYInd - rangeY - 1 + deltaY;
        dpX(pt) = dpXInd - rangeX - 1 + deltaX;
        
        processed(pt) = true;
        
        regionSet(pt) = region;
        
        
        %check all four points, if quality(pt) > quality(neighbor)
        %change quality/iniDp
        
        %up
        if quality(pt) > quality(pt - 1) && processed(pt - 1) == false
            tempQuality(pt - 1) = quality(pt);
            iniDpY(pt-1) = round( dpY(pt) );
            iniDpX(pt-1) = round( dpX(pt) );
            regionSet(pt - 1) = regionSet(pt);
        end
        
        
        %down
        if quality(pt) > quality(pt + 1)  && processed(pt + 1) == false
            tempQuality(pt + 1) = quality(pt);
            iniDpY(pt+1) = round( dpY(pt) );
            iniDpX(pt+1) = round( dpX(pt) );
            regionSet(pt + 1) = regionSet(pt);
            
        end
        
        
        %left
        if quality(pt) > quality(pt -numY) && processed(pt - numY) == false
            tempQuality(pt -numY) = quality(pt);
            iniDpY(pt - numY) = round( dpY(pt) );
            iniDpX(pt - numY) = round( dpX(pt) );
            regionSet(pt - numY) = regionSet(pt);
        end
        
        
        %right
        if quality(pt) > quality(pt + numY) && processed(pt + numY) == false
            tempQuality(pt + numY) = quality(pt);
            iniDpY(pt + numY) = round( dpY(pt) );
            iniDpX(pt + numY) = round(dpX(pt) );
            regionSet(pt + numY) = regionSet(pt);
        end
        
        region = region + 1;
        ind = ind + 1;
        
    end
end





%PERFORM BLOCK MATCHING FOR ALL OTHER POINTS
for dummyInd = totSeeds:totPoints
    
    
    [~, pt] = max(tempQuality(:) );
    [yy, xx] = ind2sub([numY, numX], pt);
    
    
    CC = 5 + normxcorr2_mex...
        (pre( kernCenterY(yy) - halfY : kernCenterY(yy) + halfY, kernCenterX(xx) - halfX: kernCenterX(xx) + halfX )...
        , post(kernCenterY(yy) - halfY - smallRangeY + iniDpY(pt) : kernCenterY(yy)  + halfY + smallRangeY + iniDpY(pt), kernCenterX(xx) - halfX - smallRangeX + iniDpX(pt): kernCenterX(xx) + halfX + smallRangeX + iniDpX(pt) ), 'valid' );
    
    
    [quality(pt), tempInd] = max(CC(:));
    tempQuality(pt) = -7;
    
    [dpYInd, dpXInd] = ind2sub(size(CC), tempInd);
    
    if dpYInd > 1 && dpYInd < 2*smallRangeY + 1
        deltaY = subSampleFit(CC(dpYInd - 1:dpYInd + 1, dpXInd) );
    else
        deltaY = 0;
    end
    
    if dpXInd > 1 && dpXInd < 2*smallRangeX + 1
        deltaX = subSampleFit(CC(dpYInd, dpXInd - 1: dpXInd + 1) );
    else
        deltaX = 0;
    end
    
    
    
    
    
    dpY(pt) = dpYInd - smallRangeY - 1 + deltaY + iniDpY(pt);
    dpX(pt) = dpXInd - smallRangeX - 1 + deltaX + iniDpX(pt);
    
    
    
    
    %%%%limit on max displacement
    if dpY(pt)  > rangeY
        dpY(pt) = rangeY;
    elseif dpY(pt) <  - rangeY
        dpY(pt) = -rangeY;
    end
    
    if dpX(pt)  > rangeX
        dpX(pt) = rangeX;
    elseif dpX(pt) <  - rangeX
        dpX(pt) = -rangeX;
    end
    %%%%%
    
    processed(pt) = true;

    
    
    
    
    
    if yy > 1
        %up
        if quality(pt) > tempQuality(pt - 1) && processed(pt - 1) == 0
            tempQuality(pt - 1) = quality(pt);
            iniDpY(pt-1) = round( dpY(pt) );
            iniDpX(pt-1) = round( dpX(pt) );
            regionSet(pt - 1) = regionSet(pt);
        end
        
    end
    
    
    if yy < numY
        
        %down
        if quality(pt) > tempQuality(pt + 1)  && processed(pt + 1) == 0
            tempQuality(pt + 1) = quality(pt);
            iniDpY(pt+1) = round( dpY(pt) );
            iniDpX(pt+1) = round( dpX(pt) );
            regionSet(pt + 1) = regionSet(pt);
            
        end
        
    end
    
    if xx > 1
        %left
        if quality(pt) > tempQuality(pt -numY) && processed(pt - numY) == 0
            tempQuality(pt -numY) = quality(pt);
            iniDpY(pt - numY) = round( dpY(pt) );
            iniDpX(pt - numY) = round( dpX(pt) );
            regionSet(pt - numY) = regionSet(pt);
        end
        
    end
    
    if xx < numX
        %right
        if quality(pt) > tempQuality(pt + numY) && processed(pt + numY) == 0
            tempQuality(pt + numY) = quality(pt);
            iniDpY(pt + numY) = round( dpY(pt) );
            iniDpX(pt + numY) = round(dpX(pt) );
            regionSet(pt + numY) = regionSet(pt);
        end
        
    end
    
    
    
end





%DROP OUT CORRECTION
processed(:) = 0;
goodPoints = ones(size(regionSet));

%%%now identify bad seeds
for r = 1:length(regionSet)
    
    
    if nnz(regionSet(:) == r ) < thresh
        
        goodPoints(regionSet == r) = 0;
        quality(regionSet == r ) = -7;
        tempQuality(regionSet == r ) = -7;
    end
end


totSeeds = 0;

for r = 1:length(seedPoints)
    pt = seedPoints(r);
    
    if goodPoints(pt) == 1;
        processed(pt) = 1;
        totSeeds = totSeeds + 1;
        %check all four points, add them back in to point set
        %up
        if  goodPoints(pt - 1) == 1 && processed(pt - 1) == 0
            tempQuality(pt - 1) = quality(pt - 1);
        end
        
        
        %down
        if goodPoints(pt + 1) == 1 && processed(pt + 1) == 0
            tempQuality(pt + 1) = quality(pt + 1);
        end
        
        
        %left
        if  goodPoints(pt - numY) == 1 && processed(pt - numY) == 0
            tempQuality(pt - numY) = quality(pt - numY);
        end
        
        
        %right
        if goodPoints(pt + numY) == 1 && processed(pt + numY) == 0
            tempQuality(pt + numY) = quality(pt + numY);
        end
        
        
    end
    
    
    
end







%now rerun the algorithm, but only change displacements at
%regions grown by bad seeds

for ind = totSeeds:totPoints
    
    
    [~, pt] = max(tempQuality(:) );
    [yy,xx] = ind2sub( [numY, numX]  , pt);
    
    if goodPoints(pt) == 0
        
        
        CC = 5 + normxcorr2_mex...
            (pre( kernCenterY(yy) - halfY : kernCenterY(yy) + halfY, kernCenterX(xx) - halfX: kernCenterX(xx) + halfX )...
            , post(kernCenterY(yy) - halfY - smallRangeY + iniDpY(pt) : kernCenterY(yy)  + halfY + smallRangeY + iniDpY(pt), kernCenterX(xx) - halfX - smallRangeX + iniDpX(pt): kernCenterX(xx) + halfX + smallRangeX + iniDpX(pt) ), 'valid' );
        
        
        [quality(pt), tempInd] = max(CC(:));
        
        [dpYInd, dpXInd] = ind2sub(size(CC), tempInd);
        
        if dpYInd > 1 && dpYInd < 2*smallRangeY + 1
            deltaY = subSampleFit(CC(dpYInd - 1:dpYInd + 1, dpXInd) );
        end
        
        if dpXInd > 1 && dpXInd < 2*smallRangeX  + 1
            deltaX = subSampleFit(CC(dpYInd, dpXInd - 1: dpXInd + 1) );
        end
        
        
        dpY(pt) = dpYInd - smallRangeY - 1 + iniDpY(pt) + deltaY;
        dpX(pt) = dpXInd - 1 - 1 + iniDpX(pt) + deltaX;
        
        
        
        
        
        
        %%%%limit on max displacement
        if dpY(pt)  > rangeY
            dpY(pt) = rangeY;
        elseif dpY(pt) <  - rangeY
            dpY(pt) = -rangeY;
        end
        
        if dpX(pt)  > rangeX
            dpX(pt) = rangeX;
        elseif dpX(pt) <  - rangeX
            dpX(pt) = -rangeX;
        end
        %%%%
        
    end
    
    tempQuality(pt) = -7;
    processed(pt) = 1;
    
    
    
    
    
    
    
    %up
    if( yy > 1 )
        
        if quality(pt) > tempQuality(pt - 1) && processed(pt - 1) == 0 && goodPoints(pt - 1) == 0
            tempQuality(pt - 1) = quality(pt);
            iniDpY(pt-1) = round( dpY(pt) );
            iniDpX(pt-1) = round( dpX(pt) );
            regionSet(pt - 1) = regionSet(pt);
            
        elseif goodPoints(pt - 1) == 1 && processed(pt - 1) == 0
            tempQuality(pt - 1) = quality(pt-1);
        end
        
    end
    
    
    %down
    if( yy < numY)
        if quality(pt) > tempQuality(pt + 1)  && processed(pt + 1) == 0 && goodPoints(pt+1) == 0
            tempQuality(pt + 1) = quality(pt);
            iniDpY(pt+1) = round( dpY(pt) );
            iniDpX(pt+1) = round( dpX(pt) );
            regionSet(pt + 1) = regionSet(pt);
            
        elseif goodPoints(pt + 1) == 1 && processed(pt + 1) == 0
            tempQuality(pt + 1) = quality(pt + 1);
        end
        
    end
    
    %left
    if( xx > 1)
        if quality(pt) > tempQuality(pt -numY) && processed(pt - numY) == 0 && goodPoints(pt - numY) == 0
            tempQuality(pt -numY) = quality(pt);
            iniDpY(pt - numY) = round( dpY(pt) );
            iniDpX(pt - numY) = round( dpX(pt) );
            regionSet(pt - numY) = regionSet(pt);
        elseif goodPoints(pt -numY) == 1 && processed(pt - numY) == 0
            tempQuality(pt - numY) = quality(pt - numY);
        end
        
    end
    %right
    if( xx < numX)
        if quality(pt) > tempQuality(pt + numY) && processed(pt + numY) == 0 && goodPoints(pt + numY) == 0
            tempQuality(pt + numY) = quality(pt);
            iniDpY(pt + numY) = round( dpY(pt) );
            iniDpX(pt + numY) = round(dpX(pt) );
            regionSet(pt + numY) = regionSet(pt);
        elseif goodPoints(pt + numY) == 1 && processed(pt + numY) == 0
            tempQuality(pt + numY) = quality(pt + numY);
        end
        
    end
    

end




quality = quality - 5;
geom = struct('startY', startY, 'stopY', stopY, 'startX', startX, 'stopX', stopX, 'stepY', stepY);




function delta=subSampleFit(vector)
delta=0;
if vector(1) ~= vector(3)
    temp= (vector(1)-vector(3))/(2*(vector(1)+vector(3)-2*vector(2)));
    if abs(temp)<1
        delta=temp;
    end
end