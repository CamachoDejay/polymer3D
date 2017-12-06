% example
clc
clear

% imaging parameter
PixelSize = 160; % nanometer

% The coefficents of the best-fit 4th order polynomial function ---------
% This parameter needs to be adjusted according to the optical imaging 
% system of the user. 
CCP = [-302.8,1386.5,-2305.3,2278.0,-1048.1]; 

% read raw image 
filename = 'bead.tif';
imRAW = imread(filename);

% perform denoising 
imDN = DeNoise(imRAW);

% perform ROI extraction
RegR = 7;      % sub-region radius 
threshold = 100;  %
[ROIs, ROI_coor] = extractROI(imRAW,imDN,RegR,threshold);

% perform 3D localization
GraR = 3; % The radius of Gradient used for caculation
spotNUM = size(ROIs,3);
spot3D = zeros(3,1);
old = load('originalOUT.mat');

for N=1:spotNUM
    [x,y,e] = GradientFit(ROIs(:,:,N),GraR);
    % reporting differences
    
    diffE = abs(e-old.e);
    fprintf('original e: %.4g ;\t new e: %.4g ;\t difference: %.4g\n',old.e,e,diffE)

    diffX = abs(x-old.x);
    fprintf('original x: %.4g ;\t new x: %.4g ;\t difference: %.4g\n',old.x,x,diffX)

    diffY = abs(y-old.y);
    fprintf('original y: %.4g ;\t new y: %.4g ;\t difference: %.4g\n',old.y,y,diffY)

    
    xc = (ROI_coor(2,N) + x)*PixelSize;
    yc = (ROI_coor(1,N) - y)*PixelSize;
    zc = CCP(1)*e^4 + CCP(2)*e^3 + CCP(3)*e^2 + CCP(4)*e + CCP(5);
    spot3D(:,N) = [xc,yc,zc];
end

% Save the localization results, [x,y,z], nanometer
save([filename(1:end-4) '_GradientFit.mat'], 'spot3D');


