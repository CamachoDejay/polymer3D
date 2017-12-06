%% 
% GradientFit.m
% Used for 3D Emitter localization by Gradient direction Fitting
% This program is free software without any warranty. This set of programs 
% is distributed as part of academic resource sharing for non-profit 
% research only. Plaese see the GNU General Public License 
% for more details.
% 
% Current implementation by Rafael Camacho (github: camachodejay). The
% following code is based on the work of Hongqiang Ma (2015) which can be
% found in http://www.pitt.edu/~liuy. Moreover the code has been heavily
% modified by refining the math and coming with cleaner definitions of
% e, x_c and y_c.
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Inputs
%   ROI  : Candidate emitter sub-region for position caculation
%          The center pixel defaults to the pixel with the peak intensity
%   RegR : Region radius of the ROI, so the size of ROI is (2*WinR)+1
%   GraR : The radius of Gradient used for caculation
%          Ignore the effect of the outer gradients with low SNR can
%          improve the performance of this method.        
% -------------------------------------------------------------------------
% Outputs
%   x,y  : The caculated emitter position in x-y dimensions
%          Unit: pixels
%   e    : The caculated ellipticity of the emitter's intensity distribution
% -------------------------------------------------------------------------
%%

function [x,y,e] = GradientFit(ROI,RegR,GraR)

% estimate the initial value of [x,y,e] by using centroid based method
[x0, y0, e0] = Centroid(ROI);

% define the coordinates of the gradient gride, set the center pixel as the original point
[m,n] = meshgrid(0.5-GraR:GraR-0.5,GraR-0.5:-1:0.5-GraR);

% define the exact gradient at each position
% I believe here they are asuming a 2D elliptical Gaussian
Gx = e0*(x0-m);
Gy = (y0+n);
G2 = (Gx.^2 + Gy.^2);
% so far G2 is a factor that only considers the ellipticity of the psf

% caculate the measured gradients - This is equation 3 in the paper
xID = RegR-GraR:RegR+GraR-1;
yID = xID;
gx = ROI(yID,xID+3)+ 2*ROI(yID+1,xID+3)+ 2*ROI(yID+2,xID+3)+ ROI(yID+3,xID+3)...
    -ROI(yID,xID)  - 2*ROI(yID+1,xID)  - 2*ROI(yID+2,xID)  - ROI(yID+3,xID);

gy = ROI(yID,xID)  + 2*ROI(yID,xID+1)  + 2*ROI(yID,xID+2)  + ROI(yID,xID+3)...
    -ROI(yID+3,xID)- 2*ROI(yID+3,xID+1)- 2*ROI(yID+3,xID+2)- ROI(yID+3,xID+3);

gx2 = gx.^2;
gy2 = gy.^2;
gxy = gx.*gy;
% g2 = (gx2 + gy2);

% Wg = sqrt(G2).*(g2);  % Wg is the weight parameter
% P = (G2.*g2)./Wg;
P = sqrt(G2);

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

A1 = (S_gxym*S_gxy - S_gy2m*S_gx2) / (S_gxy^2 - S_gy2*S_gx2);
B1 = (S_gxyn*S_gx2 - S_gx2n*S_gxy) / (S_gxy^2 - S_gy2*S_gx2);

A = (1/S_gxy)*(S_gy2m2*S_gxy - S_gxym*S_gy2m - A1*(S_gy2m*S_gxy - S_gxym*S_gy2));

Bzero = 0; % !!!!!! I can prove that B should be 0

C = (-S_gxymn) + S_gxyn*(S_gxym/S_gxy) - B1*(2*S_gy2m - S_gy2*(S_gxym/S_gxy) - S_gy2*A1);

%--------------------------------------------------------------------------

% e looks like the quadratic formula for ax^2 + bx + c = 0;
% x = [-b (+-) sqrt(b^2 - 4*a*c)]/2a
% where C = b, A = a and B = c

% the solution presented in the original code only works if C is negative.
% If not e will be 0;
errMsg = ['problems with C, must be negative.'...
          ' I have to check this due to the way the program'...
          ' was writen initially'];
assert(C<0,errMsg)
% e = (-C+sqrt(C^2-4*A*Bzero)) / (2*A);
e = -C/A;

x = A1+B1/e;

% y = A2*e+B2;
y = (-1/S_gxy) * ((-S_gxyn) - S_gy2*B1 +e*(S_gy2m -S_gy2*A1));

% Now as I did not like the way the code was writen I decided to solve the
% math on my side and come with my own solutions for the system of
% equations. These solutions are bellow and are compared with that obtained
% by the original code.

% first we calculate e
eNum = S_gy2*(S_gxym*S_gx2n - S_gxymn*S_gx2)+...
       S_gy2m*(S_gxyn*S_gx2 - S_gxy*S_gx2n) +...
       S_gxy*(S_gxymn*S_gxy - S_gxyn*S_gxym);
   
eDen = S_gy2*(S_gxym^2 - S_gy2m2*S_gx2)+...
       S_gy2m*(S_gy2m*S_gx2 - S_gxy*S_gxym) +...
       S_gxy*(S_gy2m2*S_gxy - S_gy2m*S_gxym);
   
e_new = eNum/eDen;

% then we calculate xc and yc
comD = (S_gy2*S_gx2-S_gxy^2);
K = (S_gx2n*S_gxy - S_gxyn*S_gx2)/comD;
L = (S_gy2m*S_gx2 - S_gxym*S_gxy)/comD;
M = (S_gx2n*S_gy2 - S_gxyn*S_gxy)/comD;
N = (S_gxy*S_gy2m - S_gxym*S_gy2)/comD;

x_new = K/e_new + L;
y_new = M+e_new * N;

%%
% now I report the differences:
diffE = abs(e-e_new);
fprintf('original e: %.4g ;\t new e: %.4g ;\t difference: %.4g\n',e,e_new,diffE)

diffX = abs(x-x_new);
fprintf('original x: %.4g ;\t new x: %.4g ;\t difference: %.4g\n',x,x_new,diffX)

diffY = abs(y-y_new);
fprintf('original y: %.4g ;\t new y: %.4g ;\t difference: %.4g\n',y,y_new,diffY)

