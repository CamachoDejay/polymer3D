% template for simple test of localization algorithm

clear 
% close all
clc

% simple sim of PSF in the ROI
imSize = 17;
pos_real = [2,2];

xid = 0:imSize-1;
yid = 0:imSize-1;
pix_size = 0.25;
sigX = .6;
sigY = .65;
xVal = xid.*pix_size;
yVal = yid.*pix_size;
pos_pix = (pos_real./pix_size) + 1;
sig = [sigX,sigY];
ROI = Misc.gaus2D(pos_real,sig,xVal,yVal); 

% ROI coor is always the center position
ROI_coor = [median(1:size(ROI,1)),median(1:size(ROI,1))];

GraR = 6; % The radius of Gradient used for caculation

% gradient fit for localization
[x,y,e,centOut] = Localization.gradFit(ROI,GraR);
% shifting system of reference back to corner pixel
gfitX = ROI_coor(1)+x;
gfitY = ROI_coor(2)+y;
cfitX = ROI_coor(1)+centOut.x;
cfitY = ROI_coor(2)+centOut.y;

% printing results
fprintf('Simulated emitter position [pix]: \t\t%.4f, \t%.4f \n', pos_pix(1), pos_pix(2))
fprintf('Gradient fitted emitter position [pix]: \t%.4f, \t%.4f \n', gfitX, gfitY)
fprintf('Centroid fitted emitter position [pix]: \t%.4f, \t%.4f \n', cfitX, cfitY)
fprintf('Simulated PSF Ellipticity: \t\t%.4f\n', sigY/sigX)
fprintf('Gradient fitted PSF Ellipticity: \t%.4f\n', e)
fprintf('Centroid fitted PSF Ellipticity: \t%.4f\n', centOut.e)

% making figure
figure(1)
subplot(1,2,1)
    surf(ROI)
    hold on
    scatter3(pos_pix(1),pos_pix(2),max(ROI(:)),'ro')
    hold off
    xlabel('x-pos')
    ylabel('y-pos')
    axis image
    view(2)

subplot(1,2,2)
    xSum = sum(ROI,1);
    xSum = xSum./max(xSum);
    ySum = sum(ROI,2);
    ySum = ySum./max(ySum);
    plot(1:imSize,xSum,'k')
    hold on
    plot([gfitX,gfitX],[0,1],'k')
    plot(1:imSize,ySum,'r')
    plot([gfitY,gfitY],[0,1],'r')
    hold off
    xlim([1,imSize])
a = gca;
a.FontSize = 18;
shg