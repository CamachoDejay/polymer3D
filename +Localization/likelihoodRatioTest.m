function [LRT, RMSD] = likelihoodRatioTest(im,sig,fitRes) 
im = double(im);
x = fitRes(2);
y = fitRes(1);
sizeX = size(im,1);
sizeY = size(im,2);
% number of pixels in each image
pixN = sizeX*sizeY;
% centroid of the window
xid   = 1:sizeX;
yid   = 1:sizeY;
cVal  = [mean(xid), mean(yid)];
% get likelihood ratio test
                    %   get gaussian model
                    %       position
                    pos = [cVal(1)+x,cVal(2)+y];
                    %       function
                    [G] = SimpleFitting.gauss2D(pos, sig, xid,yid,1);
                    G = G';
                    
                    %   calculation of G_bar and sGbar2.
                    %       Mean of the gaussian PSF
                    G_mean  = mean(G(:));
                    %       G bar: gaussian PSF minus its mean
                    G_bar = G - G_mean;
                    %       sum of G bar squared
                    sGbar2 = sum(sum(G_bar.^2));
                    %       convolution of the image with G_bar
                    imConv = sum(sum(G_bar.*im));
                    %       calculation of I_hat
                    I_hat  = imConv/sGbar2;
                    %       denominator, which depends on how good it fits
                    %       to random noise
                    den = var(im(:))*pixN;
                    %   calculation of the difference betwen the likelihood
                    %   of having the data explained by random noise or by
                    %   a gaussian. L(H_0)-L(H_1)
                    LRT = ( (pixN)/2 ) * (log(1-((I_hat^2 * sGbar2)/den)));
                   
                    
                    % get fit to gaussian
                    bg = ones(size(G));
                    % where A are the separate spectrum and b is what what we measured
                    meas   = im(:);
%                     meas   = meas - mean(meas);
%                     meas   = meas ./ std(meas);
                    shapes = [G(:), bg(:)];
                    
                    % find the coefficients that solve the system of equations
                    % [x, resN] = lsqnonneg(shapes,meas);
                    eps = mldivide(shapes,meas);
                    fit = shapes*eps;
                    % fit = reshape(fit,size(im));
                    dif = meas-fit;
                    RMSD = sqrt( sum(dif.^2)/length(dif));
                    