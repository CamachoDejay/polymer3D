function [x,y,e,centOut] = gradFit(ROI,GraR)
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

% Changes in calculation: Gradient Calculated twice.
%   After testing the code using PSF simulations (see test2DGradTracking.m)
%   we realized that the ellipticity is defined assymetrically. As a
%   consequence, X and Y positions which are determined from ellipticity
%   have assymetrical errors as well.
%   To account for this we now calculate the gradient twice:
%       #1 We determined X and the ellipticity
%       #2 We determined Y and the ellipticity
%   A final step will choose the ellipticity value used using 3 rules:
%       a) If one ellipticity is above 1 and the other below, we take the
%       smallest one (=elongation along X ==> more accurate)
%       b) If both are below 1 we take the highest value (closest to one)
%       c) If both are above 1 we take the smallest value (closest to one)

assert(isa(ROI,'double'),'ROI image must be a double')
% TODO: should the code not work on ints instead?

roiSize = size(ROI);
assert(all(mod(roiSize,2)==[1,1]),'input ROI dimentions must be odd')

% estimate the initial value of [x,y,e] by using centroid based method
[centOut.x, centOut.y, centOut.e] = centLocalize(ROI);

elip = zeros(2,1);
xLoc = zeros(2,1);
yLoc = zeros(2,1);

% define the coordinates of the gradient grid, set the center pixel as the original point
[m,n] = meshgrid(-GraR:GraR,-GraR:1:GraR);
%[m,n] = meshgrid(-GraR:GraR,GraR:-1:-GraR);

% caculate the measured gradients - This is equation 3 in the paper
diaShift = -2;%shift use for indices
roiXc = median(1:size(ROI,1));
roiYc = median(1:size(ROI,2));
xID = roiXc-GraR:roiXc+GraR;
yID = roiYc-GraR:roiYc+GraR;
xID = xID+diaShift;
yID = yID+diaShift;
assert(max(xID+4)<=size(ROI,1),'decrease gradient radius')
assert(max(yID+4)<=size(ROI,2),'decrease gradient radius')

gx = ROI(yID,xID+4)+ 2*ROI(yID+1,xID+4)+ 3*ROI(yID+2,xID+4)+ 2*ROI(yID+3,xID+4)+ ROI(yID+4,xID+4)...
    -ROI(yID,xID)  - 2*ROI(yID+1,xID)  - 3*ROI(yID+2,xID)  - 2*ROI(yID+3,xID) - ROI(yID+4,xID);

gy = ROI(yID,xID)  + 2*ROI(yID,xID+1)  + 3*ROI(yID,xID+2)  + 2*ROI(yID,xID+3) + ROI(yID,xID+4)...
    -ROI(yID+4,xID)- 2*ROI(yID+4,xID+1)- 3*ROI(yID+4,xID+2)- 2*ROI(yID+4,xID+3)- ROI(yID+4,xID+4);
% I changed the way n is defined, but to keep signs correct to previous
% code from Hongqiang I have to change the sign of gy here
gy = -gy;

for i = 1:2
    
    switch i
        case 1 % we do not rotate the ROI
            % guess values of x y and elip to calculate weight parameter.
            x0 = centOut.x;
            y0 = centOut.y;
            e0 = centOut.e;
            gxTmp = gx;
            gyTmp = gy;
            
        case 2 % we do ROTATE the ROI
            % guess values of x y and elip to calculate weight parameter.
            x0 = centOut.y;
            y0 = centOut.x;
            e0 = 1/centOut.e;
            gxTmp = gy';
            gyTmp = gx';
                  
        otherwise
            error('unexpected behaviour')
    end
    
    % calculating the weighting factor, in fact 1/Wg
    %%%%% USER INPUT %%%%%
    wUs = true;
    %%%%%%%%%%%%%%%%%%%%%%
     
    % define some sort of epsilon corrected eucledian distance, see coments
    % below. These numbers are only used to calculate the weighting factor,
    % Gx and Gy are only used to calculate G2.
    Gx = e0*(x0-m);
    Gy = (y0-n);
    G2 = (Gx.^2 + Gy.^2);
         
    % Now P only depends on x0 y0 and e0 so we should be able to alculate
    % it outside of the loop. However, it not simetric after doing
    % ROI'because of the way Gx is calculated. imagine that e is smaller
    % than 0, then the first time Gx will a matrix multiplied by a number
    % <0 while the second time it will be a matrix multiplied by a number
    % >0, this means that the sum(P(:)) will be different for each loop
    % iteration. This should change the overal numbers but not the
    % calculation done after as they are ratios. 
    
    if wUs
        euDist = sqrt((x0-m).^2 + (y0-n).^2);
%         g2 = (gx2 + gy2);
%         Wg = euDist.*(g2);  % Wg is the weight parameter
%         if any(Wg(:)==0)
%             warning('There is a 0 Wg value, x0 and/or y0 are integers')
%             Wg(Wg==0) = 1e-10;
%         end
%         
%         P = (G2.*g2)./Wg;
        % as g2 cancels out is better to calculate this way. NOTE that
        % sqrt(G2) is very close to an eucledian distance corrected by the
        % epsilon value
        P = (G2)./euDist; 
    else
        % There is a problem with their definition of W, in the paper Wg
        % has no e0 in the sqrt(). What is correct paper or code? probably
        % best to mail them.
        
        % NOTE that in this simple deffinition and for e0=1, P is nothing
        % else than the eucledian distance between the ROI pixel and the
        % guess for the center of the molecule [x0, y0]
        P = sqrt(G2);
    end
    
    if any(P(:)==0)
        warning('There is a 0 P value, substituting by 1e-10')
        P(P==0) = 1e-10;
    end
    
    
    % define some usefull numbers used below, so eq are easier to read
    gx2 = gxTmp.^2;
    gy2 = gyTmp.^2;
    gxy = gxTmp.*gyTmp;
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
    
    % to understand the equations go to:
    % https://camachodejay.github.io/GradientFit/
    
    % first we calculate elipticity (el)
    eNum = S_gy2*(S_gxym*S_gx2n - S_gxymn*S_gx2)+...
        S_gy2m*(S_gxyn*S_gx2 - S_gxy*S_gx2n) +...
        S_gxy*(S_gxymn*S_gxy - S_gxyn*S_gxym);
    
    eDen = S_gy2*(S_gxym^2 - S_gy2m2*S_gx2)+...
        S_gy2m*(S_gy2m*S_gx2 - S_gxy*S_gxym) +...
        S_gxy*(S_gy2m2*S_gxy - S_gy2m*S_gxym);
    
    el = eNum/eDen;
    
    % then we calculate xc and yc
    comD = (S_gy2*S_gx2-S_gxy^2);
    K = (S_gx2n*S_gxy - S_gxyn*S_gx2);
    L = (S_gy2m*S_gx2 - S_gxym*S_gxy);
    M = (S_gx2n*S_gy2 - S_gxyn*S_gxy);
    N = (S_gxy*S_gy2m - S_gxym*S_gy2);

    xLoc(i) = (K/el + L)/comD;
    yLoc(i) = (M + el*N)/comD;
    elip(i) = el;
    
end

x = xLoc(1);
y = xLoc(2);

e = elip(1);%only survives if both are equal to one

if and(elip(1)<1,elip(2)<1)
    %warning('both epsilon are smaller than one, picking the largest');
    e = max(elip(1),elip(2));
end

if and(elip(1)>1,elip(2)>1)
    % warning('both epsilon are larger than one, picking the smallest');
    e = min(elip(1),elip(2));
end

if and(elip(1)<1,elip(2)>1)
    e = elip(1);
end

if and(elip(1)>1,elip(2)<1)
    e = 1/elip(2);
end

