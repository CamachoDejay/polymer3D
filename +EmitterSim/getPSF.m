function [ psf ] = getPSF( Xgrid, Ygrid, xpos, ypos, model )
%getPSF calculates the PSF of an emitter using the desired model
%   Detailed explanation goes here

if strcmp(model.name, 'gaussian')
    sigma_x = model.sigma_x;
    sigma_y = model.sigma_y;
    

    f1 = ((Xgrid - xpos).^2) / (2*(sigma_x.^2));
    f2 = ((Ygrid - ypos).^2) / (2*(sigma_y.^2));
    A  = 1/(2*pi*sigma_x*sigma_y);  
    psf = A * exp(-(f1+f2));

end

