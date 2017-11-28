% Estimate the initial value of [x,y,e] by using centroid based method
% -------------------------------------------------------------------------
% Inputs
%   ROI    : Candidate emitter sub-region for position estimation
% -------------------------------------------------------------------------
% Outputs
%   x0,y0  : The estimated emitter position in x-y dimensions
%          
%   e0     : The estimated ellipticity of the emitter's intensity distribution
% -------------------------------------------------------------------------

function [x0, y0, e0] = Centroid(ROI)
% number of pixels in the ROI
[Ny, Nx] = size(ROI);
% radius of the ROI
w = (Ny-1)/2;
% vectors holding pixel indices
X = 1:Nx;
Y = 1:Ny;
% integral over the fisrt dimension
sIx = sum(ROI,1);
% integral over the second dimension
sIy = sum(ROI,2);
% remove off set so that the integral goes from 0 - val
sIx = sIx - min(sIx(:));
sIy = sIy - min(sIy(:));
% Centroid of an arnitrady 1D function
Xc = sum(sIx.*X)/sum(sIx);
Yc = sum(sIy.*Y')/sum(sIy);
% Create a vector of the pixel indices shifted acording to the centroid,
% this way the centroid is now at (0,0)
%   This is the same as doing X-Xc
X1 = 1-Xc:Nx-Xc;
%   This is the same as doing Y-Yc
Y1 = 1-Yc:Ny-Yc;
% What is W???? it seems to be a measurement of the width of the
% distribution. It could be realted to the second moment or quadrupole
Wx = sum(sIx.*abs(X1))/sum(sIx);
Wy = sum(sIy.*abs(Y1'))/sum(sIy);

% for i=1:(2*w+1)
%     if i<Xc-3*Wx || i>Xc+3*Wx
%         sIx(i) = 0;
%     end
%     if i<Yc-3*Wy || i>Yc+3*Wy
%         sIy(i) = 0;
%     end
% end

% they seem to be filtering some data out of the ROI Ws are some sort of
% window. Here they set the tails of the sum to 0 according to the Ws
sIx(1:floor(Xc-3*Wx)) = 0;
sIx(floor(Xc+3*Wx)+1:end) = 0;

% here there seems to be an error in the code, probably here they wanted to
% also clean the sIy and made a typo
sIx(1:floor(Yc-3*Wy)) = 0;
sIx(floor(Yc+3*Wy)+1:end) = 0;

% here they seems to refine their centroid
Xc = sum(sIx.*X)/sum(sIx);
Yc = sum(sIy.*Y')/sum(sIy);
% Create a vector of the pixel indices shifted acording to the centroid,
% this way the centroid is now at (0,0)
X1 = 1-Xc:Nx-Xc;
Y1 = 1-Yc:Ny-Yc;
% here we go again with the W
Wx = sum(sIx.*abs(X1))/sum(sIx);
Wy = sum(sIy.*abs(Y1'))/sum(sIy);
% storing the centroid value asuming that the center of the ROI is 0,0
x0 = Xc-w-1;
y0 = Yc-w-1;
% e0 is the estimated ellipticity - but, what are the Ws?
e0 = Wy/Wx;

