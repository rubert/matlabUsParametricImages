function [dpY, dpX, quality, geom ] = multi_level_2D(pre, post, searchParams)
%%Pre and post are identically sized arrays for motion tracking.
%%Search Params is a structure with the following fields:
%     rangeY        => pixels, only search range at largest kernel size
%     rangeX        => pixels, only search range at largest kernel size
%     overlap       => percentage, .75 gives a 75% overlap
%     multiWindowY       => pixels, must be an odd number, to have a well
%                       defined center pixel.  Array with four entries.  Kernel size at all 4 levels.
%     multiWindowX       => pixels, must be an odd number.  Array with four entries.  Kernel size at all 4 levels.



%output:
%dpY => displacement axially in pixels
%dpX => displacement laterally in pixels/A-lines
%quality =>  value of cross-correlation
%spacing =>  number of pixels separating location of each dp estimate
%                need it to calculate strain: (dpY(n) - dpY(n-1) )/spacing
%geometry => start and stop for making plots

%first get the envelope data
preEnvelope = abs(hilbert(pre));
postEnvelope = abs(hilbert(post));



%Going to make rangeY 15% of the kernel size for all levels past first

rangeY = round(searchParams.multiWindowY*.15);
rangeY(1) = searchParams.rangeY;

if rangeY(4) < 2;
    rangeY(4) = 2;
end

rangeX = [3,3,1,1];
rangeX(1) = searchParams.rangeX;

halfY = (searchParams.multiWindowY - 1)/2;  %assumes odd number of pixels in kernel
halfX = (searchParams.multiWindowX - 1)/2;


%%%%%%%%%%boundary considerations%%%%%%%%%%%%%%%%%%%%%
%work out the boundaries of the displacement grids so tracking doesn't go
%out of bounds, there is no zero padding
rfY =  size(pre, 1);
rfX = size(pre,2);

maxRangeY = sum(searchParams.rangeY);
maxRangeX = sum(searchParams.rangeX);

stepY = round(searchParams.multiWindowY*(1-searchParams.overlap) );
startYrf = zeros(4,1);
stopYrf = zeros(4,1);
startX = zeros(4,1);
stopX = zeros(4,1);
numY = zeros(4,1);
numX = zeros(4,1);

kernCenterY = cell(4,1);
kernCenterX = cell(4,1);

for l = 1:4
    
    if l == 1
        startYrf(l) = 1 + maxRangeY + halfY(1);
        startX(l) = 1 + maxRangeX + halfX(1);
        stopYmax = rfY - halfY(1) - maxRangeY;
        
        
        numY(l) = floor( (stopYmax - startYrf(l)) / stepY(1) ) + 1;
        
        
        stopYrf(l) = startYrf(l) + stepY(l)*(numY(l)-1);
        stopX(l) = rfX - halfX(1) - maxRangeX;
        
        numX(l) = stopX(l) - startX(l) + 1;
        
        kernCenterY{l} = startYrf(l):stepY(l):stopYrf(l);
        kernCenterX{l} = startX(l):stopX(l);
        
        
    else
        
        startYrf(l) = startYrf(l-1) + 1;
        startX(l) = startX(l-1) + 1;
        
        stopYmax = kernCenterY{l-1}(end) - halfY(l);
        numY(l) = floor( (stopYmax - startYrf(l)) / stepY(l) ) + 1;
        
        stopYrf(l) = startYrf(l) + stepY(l)*(numY(l)-1);
        stopX(l) = rfX - halfX(1) - maxRangeX;
        
        numX(l) = stopX(l) - startX(l) + 1;
        
        kernCenterY{l} = startYrf(l):stepY(l):stopYrf(l);
        kernCenterX{l} = startX(l):stopX(l);
        
    end
    
end







geomRf = struct( 'stopY', stopYrf, 'startY', startYrf, 'stepY', stepY, 'startX', startX, 'stopX', stopX );









%now I know my geometry
%perform cross correlations
%correct errors
%interpolate initial displacements for the next level, loop and repeat

for l = 1:4
    
    if l == 1
                
        iniDpY = zeros(numY(1), numX(1));
        iniDpX = zeros(numY(1), numX(1) );
       
    end
        dpY = zeros(numY(l) , numX(l) );
        dpX = zeros(numY(l) , numX(l) );
        quality = zeros(numY(l), numX(l) );
    
    for y = 1:numY(l)  %y,x are coordinates on dp grid
        for x = 1:numX(l)
            
           
            
               %use envelope data with largest kernel size,
               %papers show envelope performs better with large kernel size
               %RF data performs better with smaller kernel size
                
               if l == 1
                   CC = normxcorr2_mex...
                       (preEnvelope( kernCenterY{l}(y) - halfY(l) : kernCenterY{l}(y) + halfY(l), kernCenterX{l}(x) - halfX(l) : kernCenterX{l}(x) + halfX(l)  )...
                       , postEnvelope(kernCenterY{l}(y) + iniDpY(y,x) - halfY(l) - rangeY(l) : kernCenterY{l}(y) + iniDpY(y,x) + halfY(l) + rangeY(l), kernCenterX{l}(x) + iniDpX(y,x) - halfX(l) - rangeX(l) : kernCenterX{l}(x) + iniDpX(y,x) + halfX(l) + rangeX(l) ), 'valid' );
                   
               else
                   CC = normxcorr2_mex...
                       (pre( kernCenterY{l}(y) - halfY(l) : kernCenterY{l}(y) + halfY(l), kernCenterX{l}(x) - halfX(l) : kernCenterX{l}(x) + halfX(l)  )...
                       , post(kernCenterY{l}(y) + iniDpY(y,x) - halfY(l) - rangeY(l) : kernCenterY{l}(y) + iniDpY(y,x) + halfY(l) + rangeY(l), kernCenterX{l}(x) + iniDpX(y,x) - halfX(l) - rangeX(l) : kernCenterX{l}(x) + iniDpX(y,x) + halfX(l) + rangeX(l) ), 'valid' );
                   
                   
                   
               end
            
            %find correlation function peak, then obtain subsample
            %displacement
            
          
            [tmpQual, tempIndex] = max(CC(:));
            [subY, subX] = ind2sub( size(CC), tempIndex );
            
            if subY >1 && subY < size(CC,1) && subX > 1 && subX < size(CC,2)     %In the middle of the CC Matrix
                deltaY = subSampleFit(CC(subY-1:subY+1, subX));
                deltaX = subSampleFit(CC(subY, subX - 1: subX + 1) ); 
            elseif subY >1 && subY < size(CC,1)  %y is in the middle, x at edge
                deltaX=0;
                deltaY = subSampleFit (CC(subY-1:subY+1, subX)');
            elseif subX > 1 && subX < size(CC,2)  %x is in the middle, y at edge
                deltaX = subSampleFit (CC(subY, subX-1:subX+1));
                deltaY=0;
            else
                deltaX=0;
                deltaY=0;
            end
            
            
            dpY(y,x) = subY - rangeY(l) - 1 + deltaY + iniDpY(y,x);  %convert index to dp and add in subsample fit
            dpX(y,x) = subX - rangeX(l) - 1 + deltaX + iniDpX(y,x);
            quality(y,x) = tmpQual;
            
        end
    end
    
    
    
    %now calculate initial dp for the next level, and reset the dp estimates
    %if l ~= 4  %as long as it isn't the final level
    
    
    
	
    if l ~= 4
  %  [m2,n2]=size(dpY);
  %  xxx={(1:m2),(1:n2)};
 %   dpY = csaps(xxx,dpY, .01,xxx);
 %   dpX = csaps(xxx,dpX, .01,xxx); 
    
	%calculate initial displacements for next level.
    [iniDpY, iniDpX] = makeIniDp(dpY, dpX, geomRf, l);  %iniDpY, DpX should have integer values.

    %get rid of any displacements that are too large for the current level
    iniDpY(iniDpY > sum(rangeY(1:l))) = sum(rangeY(1:l));
    iniDpY(iniDpY < -sum(rangeY(1:l))) = -sum(rangeY(1:l));
    iniDpX(iniDpX > sum(rangeX(1:l)) ) = sum(rangeX(1:l));
    iniDpX(iniDpX < -sum(rangeX(1:l)) ) = -sum(rangeX(1:l));
    
    
    
    end
  
end




spacing = stepY(4);

geom = struct('startY', geomRf.startY(4), 'stopY', geomRf.stopY(4), 'startX', geomRf.startX(4), 'stopX', geomRf.stopX(4), 'stepY', spacing);


function delta=subSampleFit(vector)
delta=0;
if vector(1) ~= vector(3)
    temp= (vector(1)-vector(3))/(2*(vector(1)+vector(3)-2*vector(2)));
    if abs(temp)<1
        delta=temp;
    end
end






function [iniDpY, iniDpX] = makeIniDp(dpY, dpX, geom, l)

%no undersampling to worry about
%geom is a structure containing information about grid start, stop, and
%spacing in terms of RF data.

%Median filter the dispalcement fields to smooth them before
% I gocalculate the initial dp field

 %Cubic spline smoothing
    [m2,n2]=size(dpY);
    xxx={(1:m2),(1:n2)};
    dpY = csaps(xxx,dpY, .01,xxx);
    dpX = csaps(xxx,dpX, .01,xxx);       


%x,y in terms of RF points on current grid.

yCurrent = geom.startY(l):geom.stepY(l):geom.stopY(l);
xCurrent = geom.startX(l):geom.stopX(l);

yNext = geom.startY(l+1):geom.stepY(l + 1):geom.stopY(l + 1);
xNext = geom.startX(l+1):geom.stopX(l+1);


%interpoalte, fill in out of bounds values with zeros
iniDpY = interp2( xCurrent, yCurrent, dpY, xNext, yNext', 'linear', mean(dpY(:)) ) ;
iniDpX = interp2( xCurrent, yCurrent, dpX, xNext, yNext', 'linear', mean(dpX(:)) ) ;



replaceY = find(isnan(iniDpY) );


for py = 1:length(replaceY)

    [yy,~] = ind2sub( size(iniDpY), py );
    good_row = isfinite(iniDpY(yy,:) );
    iniDpY(replaceY(py) ) = mean(iniDpY(yy, good_row) );

end


replaceX = find(isnan(iniDpX) );


for px = 1:length(replaceX)

    [yy,~] = ind2sub( size(iniDpX), px );
    good_row = isfinite(iniDpX(yy,:) );
    iniDpX( replaceX(px) ) = mean(iniDpX(yy, good_row) );

end

iniDpY = round(iniDpY);
iniDpX = round(iniDpX);

