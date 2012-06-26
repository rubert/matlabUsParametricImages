function [beta, refSpecs, sampleSpecs, centersY, centersX, attenCenterBlocksY, specFreq] =...
    calculateAttenuationRPM...
    (sampleRF, refRF, lowFreq, highFreq, samplingFreq, betaRef, blockSizeMM,...
      overlap, attenuationKernelMM)

%{
    Description:
      Given a sample and reference RF data set calculate the
      attenuation coefficient using the reference phantom method.
      The output of this code is mostly used to be passed on to
      a method to calculate scatterer sizes.
      
      
    Example Usage:
    [refSpecs, sampleSpecs, centersY, centersX, attenCenterBlocksY]...
    = calculateAttenuationRPM...
    (A0, A0sample, 1, 9, 60, .5, [8,8],...
      [.85,.85],14);
      
      
    Input:
    sampleRF:  The RF data for the sample
    refRF:  The RF data of the reference phantom.  Stored in a cell
      array.  Each entry in the array is a separate image.
    lowFreq:  (MHz) The low frequency to look at when computing the czt
    highFreq:  (MHz)  The high frequency to look at when computing the czt
    samplingFreq:  (MHz)  The sampling frequency of the RF data
    betaRef:  (dB/cm MHz)  The attenuation slope of the reference phantom
    blockSizeMM:  (mm), a two element array containing the block size
    of the sample in format [axial, lateral]
    overlap:  (0-1), a two element array containing the overlap between adjacent
      blocks in the format [axial, lateral]
    attenuationKernelMM: (mm)  The least squares line fit kernel in mm

    
    output:
    attenuation:  The attenuation slope of the sample in (dB/cm MHz)
                    at a lower resolution than the RF data
      
    refSpecs:  The reference spectra.  This is an output to be used
      in spectrum calculations for 
      
    sampleSpecs:  The sample spectra.  This is an output for use
      in scatterer size estimation
      
     attenCenterPixelsY:  The RF pixels corresponding to
      the attenuation kernel center
      
     attenCenterPixelsX:  The RF A-lines correpsonding to
      the attenuation kernel centers
      
      attenCentersBlocksY:  The block centers corresponding
      to each attenuation kernel center.  So 
      sampleSpecs{ attenCenterBlocksY(y), x}
      will give the spectrum centered at attenuation(y,x)
      
      specFreq:

      
%}

      
deltaY = 1540/ (2*samplingFreq*10^6)*10^3;
deltaX = .2;
[rfPoints, rfLines] = size(sampleRF);
numRefImages = length(refRF);

%block size and position parameters
blockPixelsY = (blockSizeMM(1)/deltaY);
blockPixelsY = blockPixelsY - mod(blockPixelsY,4);
blockStepY = round(blockPixelsY*(1-overlap(1)) );
deltaZBlocksCm = blockStepY*deltaY/10;

halfBlockY = blockPixelsY/2;
startY = 1 + halfBlockY;
centersY = startY:blockStepY:rfPoints - halfBlockY;

blockPixelsX = round(blockSizeMM(2)/deltaX);
if ~mod(blockPixelsX,2)
    blockPixelsX = blockPixelsX + 1;
end

blockStepX = round(blockPixelsX*(1-overlap(2) ) );
halfBlockX = floor(blockPixelsX/2);

startX = 1 + halfBlockX;
centersX = startX:blockStepX:rfLines-halfBlockX;
%czt parameters
winLength = blockPixelsY/2;
deltaF = (highFreq - lowFreq)/winLength;
specFreq = (0:blockPixelsY/2 -1)*deltaF + lowFreq;
percentUnitCircle = (highFreq - lowFreq)/samplingFreq;

w = exp( -(1j*2*pi/winLength)*percentUnitCircle);
a = exp( (1j*2*pi)*lowFreq/samplingFreq);

%frequency smoothing parameters
frequencySmoothingKernelSizeMHz = .25;
frequencySmoothingKernelSizePixels = floor(frequencySmoothingKernelSizeMHz/deltaF);
if frequencySmoothingKernelSizePixels < 2;
    frequencySmoothingKernelSizePixels = 2;
end

%attenuation kernel size
attenuationKernelBlocks = round( attenuationKernelMM/(blockStepY*deltaY) );



%%Calculate the reference spectrum at each depth and 
%%average it laterally
refSpecs = cell(length(centersY),1 );
for y  =1:length(centersY)
    refSpecs{y} = zeros(winLength,1);
end

for im = 1:numRefImages
    tempRef = refRF{im};
    for y =1:length(centersY)
            refBlock = tempRef(centersY(y)-halfBlockY:centersY(y) + halfBlockY, :);
            tempSpec = calculateSpectrum(refBlock,w,a);
            refSpecs{y} = refSpecs{y} + tempSpec;
    end
    
end

for y = 1:length(centersY)
    refSpecs{y} = refSpecs{y}/max(refSpecs{y});
end


sampleSpecs = cell(length(centersY), length(centersX));

%%%Next loop through the A-lines, calculating Attenuation slope as I go
lnRS = cell(length(centersY),1);
freqSlope = zeros(attenuationKernelBlocks,1);
beta = zeros(length(centersY) - attenuationKernelBlocks, length(centersX) );


attenCenterBlocksY = (1:(length(centersY) - attenuationKernelBlocks)) + floor(attenuationKernelBlocks/2);

for x = 1:length(centersX)
    
    %get the ln(RS(f,z))
    for y =1:length(centersY)
        sampleBlock = sampleRF(centersY(y)-halfBlockY:centersY(y) + halfBlockY, centersX(x)-halfBlockX:centersX(x)+halfBlockX);
        tempSampleSpec = calculateSpectrum(sampleBlock,w,a);
        tempSampleSpec = conv(tempSampleSpec, ones(frequencySmoothingKernelSizePixels,1), 'same')/frequencySmoothingKernelSizePixels;
        sampleSpecs{y,x} = tempSampleSpec;
        lnRS{y} = log(tempSampleSpec./refSpecs{y});
    end
    
    
    %over an attenuation kernel get regression of logRS at a frequency
    %then regression of slope with depth
    for y = 1:length(centersY) - attenuationKernelBlocks
        
        for b = 0:attenuationKernelBlocks-1
           freqSlope(b+1) = lsqFit(lnRS{y+b}, deltaF);
        end
        slope = lsqFit(freqSlope, deltaZBlocksCm);
        beta(y,x) = -8.686*slope/(4) + betaRef;
    end
    
end


function slope = lsqFit(shifts, deltaZ)
%{
  Uses the shift between blocks (Mhz) in conjuction
  with the block spacing (cm) to calculate 
  dfc/dz

%}

%%%
%matrix equation is
% 
%      [ z1   1    [ a      =  [dp1
%        z2   1      b ]        dp2
%        z3   1 ]               dp3 ]

position = deltaZ*(1:length(shifts))';
A = [position, ones(length(shifts),1)];

        
b = shifts;
temp = lscov(A,b);
slope = temp(1);
        

        
function spec = calculateSpectrum(block, w,a)
%{
    Compute a fourier spectrum by averaging over 3 depths and over A-lines.

    Input:
    block:  Block of RF data to do a Welch bartlett style estimation on.
            50% window overlap.  50% window size
    w:  CZT parameter
    a:  CZT parameter

%}

    %make axial block length divisible by 4
    [points, lines] = size(block);
    points = points - mod(points,4);
    block = block(1:points, :);
    avg = zeros(points/2,1);
    for ax = 1:3
        for l = 1:lines
            tempFFT = czt(hanning(points/2).*block(1 + (ax-1)*points/4:1 + (ax-1)*points/4 + points/2 - 1, l), points/2, w, a);
            avg = avg + abs(tempFFT).^2;
        end
    end
    
    spec = avg/max(avg(:));
