function [ BW, L0mL1, FAR ] = GLRTfiltering( imIn, delta, FWHM_pix, chi2  )
%GLRTfiltering does a Generalized Likelihood Ratio Test on the input image.
%The null hypothesis being that the pixel and its sourranding can be
%explained by noise, and the test hypothesis that it can be explained by a
%gaussian signal.
%   Based on "Dynamic multiple-target tracing to probe spatiotemporal
%   cartography of cell membranes" by Sergé et. al.
%   This code uses the vectorization capabilities of matlab, this helps
%   speed but might make the code more difficult to read. to understand the
%   notation please go to the SI of the article referred above.
% 
%   author: Rafael Camacho, camachodejay@yahoo.com

% testing for preconditions:
assert(length(size(imIn))==2, 'input image must be a m-by-n matrix')
assert(isa(imIn,'double'),'input image must be a double, if it is integers calculations will fail')
assert(delta < min(size(imIn))/2, 'your delta (window size) is too big')
assert(all([isreal(delta), rem(delta,1)==0, delta>0]),...
    'delta must be a real positive whole number eg. 3.0')
assert(chi2>0, 'chi2 used for constant false alarm rate test must be postive')

% size of the window (ROI) used for testing, if we test pixel [5,4] then we
% use delta as the radius around the pixel to be used for test. Moreover
% our test is based on a squared ROI.
Ws = 1+delta*2;

% calculations for estimating the PSF of the emitter -G-, which we
% approximate by a gaussian.
%   calculating the sigma of the gaussian from the FWHM
sigma_pix = FWHM_pix / sqrt(8*log(2));
%   calculating a gaussian PSF located in the middle of the ROI
psf_model.name = 'gaussian';
psf_model.sigma_x = sigma_pix;
psf_model.sigma_y = psf_model.sigma_x;
[Xgrid, Ygrid] = meshgrid(-delta:1:delta,-delta:1:delta);
[ G ] = emitter_psf( Xgrid, Ygrid, 0, 0, psf_model);

% calculation of G_bar and sGbar2.
%   Mean of the gaussian PSF
G_mean  = mean(G(:));
% G bar: gaussian PSF minus its mean
G_bar = G - G_mean;
% sum of G bar squared
sGbar2 = sum(sum(G_bar.^2));
% the above is equivalent to do:
%   G_mean  = mean(G(:));
%   G_sqr_mean = mean(G(:).^2);
%   G_var   = G_sqr_mean - G_mean^2;
%   sGbar2 = G_var*(Ws^2);

% size of the input image
s1 = size(imIn,1);
s2 = size(imIn,2);
% size of the usefull area
ud_s1  = s1-(2*delta);
ud_s2  = s2-(2*delta);

% calculation of sXbar2 - Sum of X bar squared - using mean filter
% mean filter of the image.
X_win_mean     = imboxfilt(imIn,Ws);
% mean filter of the image squared
X_sqr_win_mean = imboxfilt(imIn.^2,Ws);
% calculating the image of the variance
X_win_var      = X_sqr_win_mean - X_win_mean.^2;
% Image that contains the sum of the variance squared
% sXbar2 is variance of image multiplied by Ws^2, the name sXbar2 might not
% be very accurate: TODO FIX NAME
sXbar2 = X_win_var.*(Ws^2);
% removing the useless part of the image
sXbar2 = sXbar2(delta+1:1:s1-delta,delta+1:1:s2-delta);
% calculation of I hat as a convolution of the input image with G bar all
% divided by the sum of g bar squared
I_hat = conv2(imIn,G_bar,'valid')./sGbar2;
% calculation of the difference betwen the likelihood of having the data
% explained by random noise or by a gaussian. L(H_0)-L(H_1)
L0mL1   = ( (Ws^2)/2 ) * (log(1-(((I_hat.^2) .* sGbar2)./sXbar2)));
% Detection at constant false alarm rate (CFAR)
FAR = -2*(L0mL1);
BW = true([ud_s1,ud_s2]);
BW(FAR > chi2) = false;
BW = (~BW);
BW = padarray(BW,[delta delta]);
FAR = padarray(FAR,[delta delta]);
% post conditions
assert(all(size(BW) == size(imIn)),'Ups, problems! size of output image makes no sense')

end

function [ psf ] = emitter_psf( Xgrid, Ygrid, xpos, ypos, model )
%EMITTER_PSF calculates the PSF using the desired model

if strcmp(model.name, 'gaussian')
    sigma_x = model.sigma_x;
    sigma_y = model.sigma_y;
    

    X1 = ((Xgrid - xpos).^2) / (2*(sigma_x.^2));
    X2 = ((Ygrid - ypos).^2) / (2*(sigma_y.^2));
    A  = 1/(2*pi*sigma_x*sigma_y);  
    psf = A * exp(-(X1+X2));

end
end

