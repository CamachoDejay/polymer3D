function [noisyROI] = generateNoise (ROI, noiseType, noiseProp)

%Assertion
assert(ismatrix(ROI),'ROI input should be a vector or a matrix');
assert(isnumeric(ROI),'ROI input should contain numbers');
assert(ischar(noiseType),'noiseType should be a string');
assert(isstruct(noiseProp), 'noiseProp should be a structure');

  switch noiseType
        case 'Gaussian'
            %Add background to avoid negative value after adding noise           
            tmpROI = ROI +noiseProp.bkg;%to avoid negative value with Gaussian noise
            stdev  = noiseProp.maxCount/noiseProp.S2N;
            gNoise = stdev*randn(size(tmpROI))+ 0;%randn yield var = 1 multiplying
            %the noise will multiply variance the same amount to the power
            %of 2, thus we sqrt the variance demanded by the user.
            ROI =tmpROI+gNoise;%Rounded to obtain integers while keeping
            %double type for gradient function.
           
        case 'Poisson'
           tmpROI = uint16(ROI +bkg); %Cast to int for imnoise
           ROI = double(imnoise(tmpROI,'poisson'));   
   
      otherwise            
  end
   
  noisyROI = double(round(abs(ROI)));%Rounded to obtain integers while keeping
            %double type for gradient function.

end