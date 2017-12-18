function [x,y,elip,centOut] = gradFit(ROI,GraR)
% GRADFIT: estimates the position [x, y] and ellipticity [elip] of an
% emitter' PSF imaged in ROI via gradient direction fitting. This function
% can be used for 3D Emitter localization.
%   Inputs:
%       ROI: Region of interest containing the centered image of a single
%       emitter.
%       GraR: The radius of gradient used for caculation. Ignorng the
%       effect of the outer gradients with low SNR can improve the
%       performance of this method.
%   Outputs:
%       x_c,y_c: Estimated position of the emitter.
%       elip:    Estimated ellipticity of the emitters PSF.   
%
%   Current implementation by Rafael Camacho (github: camachodejay). The
%   following code is based on the work of Hongqiang Ma (2015) which can be
%   found in http://www.pitt.edu/~liuy. Moreover the code has been heavily
%   modified by refining the math and coming with cleaner definitions of
%   e, x_c and y_c.

assert(isa(ROI,'double'),'ROI image must be a double')
% TODO: should the code not work on ints instead? 

roiSize = size(ROI);
assert(all(mod(roiSize,2)==[1,1]),'input ROI dimentions must be odd')

% estimate the initial value of [x,y,e] by using centroid based method
[x0, y0, e0] = centLocalize(ROI);
centOut.x = x0;
centOut.y = y0;
centOut.e = e0;
% define the coordinates of the gradient grid, set the center pixel as the original point
[m,n] = meshgrid(-GraR:GraR,-GraR:1:GraR);

% define the exact gradient at each position. I believe here they are
% asuming a 2D elliptical Gaussian, these numbers are only used to
% calculate the weighting factor
Gx = e0*(x0-m);
Gy = (y0-n);
G2 = (Gx.^2 + Gy.^2);

% caculate the measured gradients - This is equation 3 in the paper
roiXc = median(1:size(ROI,1));
roiYc = median(1:size(ROI,2));
xID = roiXc-GraR:roiXc+GraR;
yID = roiYc-GraR:roiYc+GraR;
xID = xID-1;
yID = yID-1;

assert(max(xID+3)<=size(ROI,1),'decrease gradient radius')
assert(max(yID+3)<=size(ROI,2),'decrease gradient radius')

gx = ROI(yID,xID+3)+ 2*ROI(yID+1,xID+3)+ 2*ROI(yID+2,xID+3)+ ROI(yID+3,xID+3)...
    -ROI(yID,xID)  - 2*ROI(yID+1,xID)  - 2*ROI(yID+2,xID)  - ROI(yID+3,xID);

gy = ROI(yID,xID)  + 2*ROI(yID,xID+1)  + 2*ROI(yID,xID+2)  + ROI(yID,xID+3)...
    -ROI(yID+3,xID)- 2*ROI(yID+3,xID+1)- 2*ROI(yID+3,xID+2)- ROI(yID+3,xID+3);

% I changed the way n is defined, but to keep signs correct to previous
% code from Hongqiang I have to change the sign of gy here
gy = -gy;
% TODO: if I use ints here we can not handle negative numbers and we should. 

gx2 = gx.^2;
gy2 = gy.^2;
gxy = gx.*gy;

% calculating the weighting factor, in fact 1/Wg
% g2 = (gx2 + gy2);
% Wg = sqrt(G2).*(g2);  % Wg is the weight parameter
% P = (G2.*g2)./Wg;
P = sqrt(G2);
if any(P(:)==0)
    warning('There is a 0 P value, substituting by 1e-10')
    P(P==0) = 1e-10;
end
% There is a problem with their definition of W, in the paper Wg has no e0
% in the sqrt(). What is correct paper or code? probably best to mail them.

% solve the equation to get the best fit [x,y,e] --------------------------
S_gy2m = sum(sum(gy2.*m./P));
S_gy2 = sum(sum(gy2./P));
S_gxy = sum(sum(gxy./P));
S_gxyn = sum(sum(gxy.*n./P));
S_gxym = sum(sum(gxy.*m./P));
S_gx2 = sum(sum(gx2./P));
S_gx2n = sum(sum(gx2.*n./P));
S_gy2m2 = sum(sum(gy2.*m.^2./P));
S_gxymn = sum(sum(gxy.*m.*n./P));

% I did not like the way the code was writen I decided to solve the
% math on my side and come with my own solutions for the system of
% equations. These solutions are consistent with previous code.

% first we calculate e
eNum = S_gy2*(S_gxym*S_gx2n - S_gxymn*S_gx2)+...
       S_gy2m*(S_gxyn*S_gx2 - S_gxy*S_gx2n) +...
       S_gxy*(S_gxymn*S_gxy - S_gxyn*S_gxym);
   
eDen = S_gy2*(S_gxym^2 - S_gy2m2*S_gx2)+...
       S_gy2m*(S_gy2m*S_gx2 - S_gxy*S_gxym) +...
       S_gxy*(S_gy2m2*S_gxy - S_gy2m*S_gxym);
   
elip = eNum/eDen;

% then we calculate xc and yc
comD = (S_gy2*S_gx2-S_gxy^2);
K = (S_gx2n*S_gxy - S_gxyn*S_gx2)/comD;
L = (S_gy2m*S_gx2 - S_gxym*S_gxy)/comD;
M = (S_gx2n*S_gy2 - S_gxyn*S_gxy)/comD;
N = (S_gxy*S_gy2m - S_gxym*S_gy2)/comD;

x = K/elip + L;
y = M+elip * N;

% now I have to add for the fact that the gradient operant is even while my
% image is odd
x = x +0.5;
y = y +0.5;

