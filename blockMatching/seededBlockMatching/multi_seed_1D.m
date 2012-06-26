function [dpY, quality, geom] = multi_seed_1D(pre, post, params)
%%
%I.  DESCRIPTION
%This code is an implementation of the seeds algorithm proposed by Lujie
%Chen et al in A quality guided displacement tracking algorithm for
%ultrasonic elasticity imaging.  The quality metric is the normalized
%cross-correlation, and the sub-sample interpolation of the cross
%correlation is done through a parabolic fit.  This algorithm performs 1-D
%motion tracking axially.
%
%
%
%
%II.  OUTPUT
%dpY - the displacements in the y direction, units are pixels
%quality - The value of the normalized CC associated with each Dp estimate
%geom - a structure containing the following fields:
%        stepY -  spacing between displacement estimates (pixels)
%        startY - location of first displacemetn estimate relative to RF data
%        stopY - location of last displacemetn estimate relative to RF data
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
%         overlap :  The overlap between adjacent blocks axially, a
%         percentage expressed as a number between 0 and 1
%
%         windowY :  The block matching kernel size in the axial
%         direction, THIS NUMBER MUST BE AN ODD INTEGER
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
smallRangeY = params.smallRangeY;
[rfPoints, numLines] = size(pre);




%work out the number of displacement points that will fit
startY = 1 + rangeY + halfY  +smallRangeY ;  %first window center
stepY = round(params.windowY*(1-params.overlap));
upperBound = rfPoints - halfY - 1 - rangeY - smallRangeY ;
numY = floor( (upperBound - startY) /stepY );


%calculate window centers
kernCenterY = zeros(numY, 1);
kernCenterY(1) = startY;

for yy = 2:numY
    
    kernCenterY(yy) = kernCenterY(yy-1) + stepY;
    
end

stopY = kernCenterY(end);

totPoints = numY*numLines;

%the seed step is to chunk the image up into squares
%having a certain percent of the image area
ySeedLength = round( .1*numY );
xSeedLength = round( .1*numLines );

thresh = xSeedLength*ySeedLength;

%making sure all the initial seeds are interior points
%to simplify code
seedsY = 2:ySeedLength:numY-1;
seedsX = 2:xSeedLength:numLines-1;
totSeeds = length(seedsY)*length(seedsX);


%region set will keep track of which seed a group of points belongs
%to

regionSet = zeros(numY,numLines);

dpY = zeros(numY, numLines);
quality = -7*ones(numY, numLines);
processed = zeros(totPoints, 1);

pointSet = zeros(totPoints,1);
iniDp = zeros(totPoints,1);

region = 1;

seedPoints = zeros(length(seedsY)*length(seedsX),1 );


%PERFORM CROSS CORRELATION FOR SEEDS
ind = 1;

for yy = 1: length(seedsY)
    for xx = 1:length(seedsX)
        yInd = seedsY(yy);
        xInd = seedsX(xx);
        pt = sub2ind( [numY, numLines], yInd, xInd);
        seedPoints(ind) = pt;
        
        CC = 5 + normxcorr2_mex...
            (pre( kernCenterY(yInd) - halfY : kernCenterY(yInd) + halfY, xInd )...
            , post(kernCenterY(yInd) - halfY - rangeY : kernCenterY(yInd)  + halfY + rangeY, xInd ), 'valid' );
        
        
        [quality(pt), tempInd] = max(CC);
        
        if tempInd > 1 && tempInd < length(CC(:))
            deltaY = subSampleFit(CC(tempInd-1:tempInd+1));
        else
            deltaY = 0;
        end
        
        dpY(pt) = tempInd - rangeY - 1 + deltaY;
        
        processed(pt) = 1;
        
        regionSet(pt) = region;
        
        %check all four points, if quality(pt) > quality(neighbor)
        %change quality/iniDp
        
        %up
        if quality(pt) > quality(pt - 1) && processed(pt - 1) == 0
            pointSet(pt - 1) = 1;
            quality(pt - 1) = quality(pt);
            iniDp(pt-1) = round(dpY(pt) );
            regionSet(pt - 1) = regionSet(pt);
        end
        
        
        %down
        if quality(pt) > quality(pt + 1)  && processed(pt + 1) == 0
            pointSet(pt + 1) = 1;
            quality(pt + 1) = quality(pt);
            iniDp(pt+1) = round( dpY(pt) );
            regionSet(pt + 1) = regionSet(pt);
            
        end
        
        
        %left
        if quality(pt) > quality(pt -numY) && processed(pt - numY) == 0
            pointSet(pt - numY) = 1;
            quality(pt -numY) = quality(pt);
            iniDp(pt - numY) = round( dpY(pt) );
            regionSet(pt - numY) = regionSet(pt);
        end
        
        
        %right
        if quality(pt) > quality(pt + numY) && processed(pt + numY) == 0
            pointSet(pt + numY) = 1;
            quality(pt + numY) = quality(pt);
            iniDp(pt + numY) = round( dpY(pt) );
            regionSet(pt + numY) = regionSet(pt);
        end
        
        region = region + 1;
        ind = ind + 1;
    end
end




%PERFORM BLOCK MATCHING FOR ALL OTHER POINTS
for ind = totSeeds:totPoints
    
    
    
    [~, pt] = max(quality(:) .*pointSet);
    [ptY,x] = ind2sub( [numY, numLines]  , pt);
    
    CC = 5 + normxcorr2_mex...
        (pre( kernCenterY(ptY) - halfY   : kernCenterY(ptY) + halfY , x )...
        , post(kernCenterY(ptY) - halfY - smallRangeY + iniDp(pt)  : kernCenterY(ptY)  + halfY + smallRangeY + iniDp(pt) , x ), 'valid' );
    
    
    
    [quality(pt), tempInd] = max(CC);
    
    if tempInd > 1 && tempInd < length(CC(:))
        deltaY = subSampleFit(CC(tempInd-1:tempInd+1));
    else
        deltaY = 0;
    end
    
    dpY(pt) = tempInd - smallRangeY - 1 + iniDp(pt) + deltaY;
    
    
    if dpY(pt)  > rangeY
        dpY(pt) = rangeY;
    elseif dpY(pt) <  - rangeY
        dpY(pt) = -rangeY;
    end
    
    iniDp(pt) = round(dpY(pt));
    
    processed(pt) = 1;
    pointSet(pt) = 0;
    
    
    
    
    
    %up
    if(ptY > 1)
        if quality(pt) > quality(pt - 1) && processed(pt - 1) == 0
            pointSet(pt - 1) = 1;
            quality(pt - 1) = quality(pt);
            iniDp(pt-1) = iniDp(pt);
            regionSet(pt - 1) = regionSet(pt);
        end
    end
    
    if( ptY < numY )
        %down
        if quality(pt) > quality(pt + 1)  && processed(pt + 1) == 0
            pointSet(pt + 1) = 1;
            quality(pt + 1) = quality(pt);
            iniDp(pt+1) = iniDp(pt);
            regionSet(pt + 1) = regionSet(pt);
            
        end
        
    end
    if(x < numLines)
        %right
        if quality(pt) > quality(pt + numY) && processed(pt + numY) == 0
            pointSet(pt + numY) = 1;
            quality(pt + numY) = quality(pt);
            iniDp(pt + numY) = iniDp(pt);
            regionSet(pt + numY) = regionSet(pt);
        end
    end
    if(x > 1)
        %left
        if quality(pt) > quality(pt -numY) && processed(pt - numY) == 0
            pointSet(pt - numY) = 1;
            quality(pt -numY) = quality(pt);
            iniDp(pt - numY) = iniDp(pt);
            regionSet(pt - numY) = regionSet(pt);
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
            pointSet(pt - 1) = 1;
        end
        
        
        %down
        if goodPoints(pt + 1) == 1 && processed(pt + 1) == 0
            pointSet(pt + 1) = 1;
        end
        
        
        %left
        if  goodPoints(pt - numY) == 1 && processed(pt - numY) == 0
            pointSet(pt - numY) = 1;
        end
        
        
        %right
        if goodPoints(pt + numY) == 1 && processed(pt + numY) == 0
            pointSet(pt + numY) = 1;
        end
        
        
    end
    
    
    
end







%now rerun the algorithm, but only change displacements at
%regions grown by bad seeds

for ind = totSeeds:totPoints
    
    
    
    [~, pt] = max(quality(:) .*pointSet);
    [ptY,x] = ind2sub( [numY, numLines]  , pt);
    
    if goodPoints(pt) == 0 
        CC = 5 + normxcorr2_mex...
            (pre( kernCenterY(ptY) - halfY   : kernCenterY(ptY) + halfY , x )...
            , post(kernCenterY(ptY) - halfY - smallRangeY + iniDp(pt)  : kernCenterY(ptY)  + halfY + smallRangeY + iniDp(pt) , x ), 'valid' );
        
        
        
        [quality(pt), tempInd] = max(CC);
        
        if( (tempInd > 1) && (tempInd < length(CC) )   )
            deltaY = subSampleFit(CC(tempInd-1:tempInd+1));
        else
            deltaY = 0;
        end
        
        dpY(pt) = tempInd - smallRangeY - 1 + iniDp(pt) + deltaY;
        
        if dpY(pt)  > rangeY
            dpY(pt) = rangeY;
        elseif dpY(pt) <  - rangeY
            dpY(pt) = -rangeY;
        end
        
        iniDp(pt) = round(dpY(pt)) ;
        
    end
    
    processed(pt) = 1;
    pointSet(pt) = 0;
    
    
    
    
    %up
    if(ptY > 1 )
        
        if quality(pt) > quality(pt - 1) && processed(pt - 1) == 0 && goodPoints(pt - 1) == 0
            pointSet(pt - 1) = 1;
            quality(pt - 1) = quality(pt);
            iniDp(pt-1) = iniDp(pt);
            regionSet(pt - 1) = regionSet(pt);
            
        elseif goodPoints(pt - 1) == 1
            pointSet(pt - 1) = 1;
        end
        
    end
    
    
    %down
    if(ptY < numY)
        
        if quality(pt) > quality(pt + 1)  && processed(pt + 1) == 0 && goodPoints(pt + 1) == 0
            pointSet(pt + 1) = 1;
            quality(pt + 1) = quality(pt);
            iniDp(pt+1) = iniDp(pt);
            regionSet(pt + 1) = regionSet(pt);
            
        elseif goodPoints(pt + 1) == 1
            pointSet(pt + 1) = 1;
            
        end
        
    end
    
    %left
    if(x > 1)
        if quality(pt) > quality(pt -numY) && processed(pt - numY) == 0 && goodPoints(pt - numY) == 0
            
            pointSet(pt - numY) = 1;
            quality(pt -numY) = quality(pt);
            iniDp(pt - numY) = iniDp(pt);
            regionSet(pt - numY) = regionSet(pt);
            
            
        elseif goodPoints(pt - numY) == 1
            pointSet(pt - numY) = 1;
            
        end
    end
    %right
    if(x < numLines)
        
        if quality(pt) > quality(pt + numY) && processed(pt + numY) == 0 && goodPoints(pt + numY) == 0
            pointSet(pt + numY) = 1;
            quality(pt + numY) = quality(pt);
            iniDp(pt + numY) = iniDp(pt);
            regionSet(pt + numY) = regionSet(pt);
            
        elseif goodPoints(pt + numY) == 1
            pointSet(pt + numY) = 1;
            
            
        end
        
    end
    
    
end




quality = quality - 5;

geom = struct;
geom.stepY = stepY;
geom.startY = startY;
geom.stopY = stopY;



function delta=subSampleFit(vector)
delta=0;
if vector(1) ~= vector(3)
    temp= (vector(1)-vector(3))/(2*(vector(1)+vector(3)-2*vector(2)));
    if abs(temp)<1
        delta=temp;
    end
end

