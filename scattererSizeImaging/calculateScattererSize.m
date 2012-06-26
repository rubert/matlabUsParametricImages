function [scatSizeImage, attenuation, scatSizeCentersY, scatSizeCentersX] = calculateScattererSize...
    (sampleRF, refRF, refBSC, refBscStart, refBscStep,lowFreq, highFreq, samplingFreq, betaRef, blockSizeMM,overlap, ...
      attenuationKernelMM)

%{
      
    Description
    This code uses the reference phantom method of Yao et al to
      estimate attenuation, then it calculates a set of
      theoretical backscatter coefficients according to
      Faran scattering theory.
      
      The backscatter curve is estimated at each block, then
      the square difference of the log of the backscatter coefficients
      is computed for a range of scatterer sizes.
      
      The scatter size with the lowest error at a given depth
      is the estimated scatter size.
      
      
    Example Usage:
   [scatSizeImage, attenuation, scatSizeCentersY, scatSizeCentersX] = calculateScattererSize...
    (A0, rfRef,bscRef, startFreqBsc, deltaFBsc, 6, 11, 40, .5, [8,10], [.85,.85], ...
      14);
 
    Input:
    sampleRF:  The RF data for the sample
    refRF:  The RF data of the reference phantom.  This should probably
      be multiple planes.  Each plane of RF data is a separate entry in a 
      cell array.
    refBSC:  The backscatter coefficients of the reference phantom 
    in the frequency range of interest
    refBscStart:(MHz)  The frequency at which the backscatter coefficient starts
    refBscStep: (MHz)  The frequency interval between Bsc's
    lowFreq:  (MHz) The low frequency to look at when computing the czt
    highFreq:  (MHz)  The high frequency to look at when computing the czt
    samplingFreq:  (MHz)  The sampling frequency of the RF data
    betaRef:  (dB/cm MHz)  The attenuation slope of the reference phantom
    blockSizeMM:  (mm), a two element array containing the block size
    of the sample in format [axial, lateral]
    blockOverlap:  (0-1), a two element array containing the overlap between adjacent
      blocks in the format [axial, lateral]
    attenuationKernelMM: (mm)  The least squares line fit kernel in mm
    focalDepthMM: (MM) The focal depth of the image.  Used to fit a gaussian
      to power spectrum
    
    output:
    scatSizeImage:  A scatterer size image in microns
    attenuation:  An attenuation map in dB/cm MHz
    scatSizeCentersY:  The location of the scatterer size and 
      attenuation image pixels relative to RF data, axially.
    scatSizeCentersX:  The location of the scatterer size and 
      attenuation image pixels relative to RF data, laterally.
      The pixel in rf(scatSizeCentersY(y), scatSizeCentersX(x))
      and scatSizeImage(y,x) have the same location.
      
%}


%First work out the attenuation image, then work out scatter sizes
 [attenuation,refSpecs, sampleSpecs, centersY, centersX, attenCenterBlocksY, specFreq] = calculateAttenuationRPM...
  (sampleRF, refRF, lowFreq, highFreq, samplingFreq, betaRef, blockSizeMM,...
      overlap, attenuationKernelMM);
  
%use interpolation if necessary to compute backscattercoefficients at
%same frequencies as power spectrum
knownBscFreq = (0:length(refBSC)-1)*refBscStep + refBscStart;
matchedBsc = interp1(knownBscFreq, refBSC, specFreq);
 
%for plotting relative to RF data  
scatSizeCentersY = centersY(attenCenterBlocksY);
scatSizeCentersX = centersX;
  
%first calculate theoretical faran backscatter curves in increments
%of one micron from 10 to 150 microns
diam = .1:.1:150;
faranBsc = cell(length(diam),1);


for counter = 1:length(diam)
  faranBsc{counter} = sphere(specFreq, diam(counter),1540,5570,3374.7,1020,2540,180,50);
end
  

deltaYCm = 1540/(2*samplingFreq*10^6)*10^2;

scatSizeImage = zeros( length(attenCenterBlocksY), length(centersX) );
%%Iterate over all points
for y = 1:length(attenCenterBlocksY)
    for x = 1:length(centersX)
        
        betaDiff = mean(attenuation(1:y,x)) - betaRef;
        betaDiffNepers = betaDiff/8.686;
        blockIdx = attenCenterBlocksY(y);
        depthCm = centersY(blockIdx)*deltaYCm;
        bscSample =matchedBsc'.*sampleSpecs{blockIdx,x}./refSpecs{blockIdx}.*exp(4*betaDiffNepers.*specFreq'*depthCm);
        scatSizeImage(y,x) = findScatterSize(bscSample, faranBsc, diam);
        
        
    end
end




function bestMatchSize = findScatterSize(bscSample, faranBsc, diam)
%{

    Compute the best scatter size according to the method in
    Gerig et al and Insana et al.

%}

mmse = zeros(length(diam), 1);

for d= 1:length(diam)
   
    psi = 10*log(bscSample) - 10*log(faranBsc{d}');
    psiBar = mean(psi);
    
    temp = (psi - psiBar).^2;
    mmse(d) = mean(temp(:));
    
end

[~, minInd] = min(mmse);

bestMatchSize= diam(minInd);