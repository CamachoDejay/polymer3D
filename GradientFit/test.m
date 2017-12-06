clear 
close all
clc

pix_size = 0.25;
pos_real = [1.5,1.8];
sigX = .6;
sigY = 1.2;


im_size = 13;
xid = 0:im_size-1;
yid = 0:im_size-1;

xVal = xid.*pix_size;
yVal = yid.*pix_size;
pos_pix = (pos_real./pix_size) + 1


sig = [sigX,sigY];
ROI = gaus2D(pos_real,sig,xVal,yVal); 

% ROI coor is always the center position
ROI_coor = [median(1:size(ROI,1)),median(1:size(ROI,1))];

GraR = 4; % The radius of Gradient used for caculation

[x,y,e] = GradientFit(ROI,GraR);

% fprintf('x out: %.4g \n',x)
% fprintf('y out: %.4g \n',y)

xc = (ROI_coor(1) + x);
yc = (ROI_coor(2) + y);
fprintf('x pos [pix]: %.4g  [fit]: %.4g [real]: %.4g \t \n',xc,(xc-1)*pix_size,pos_real(1))
fprintf('y pos [pix]: %.4g  [fit]: %.4g [real]: %.4g \t \n',yc,(yc-1)*pix_size,pos_real(2))
fprintf('elipt: %.4g \n',e)


figure(1)
% surf(xid,yid,G)
subplot(1,2,1)
surf(ROI)
% contourf(xid,yid,G)
xlabel('x-pos')
ylabel('y-pos')
axis image
view(2)
subplot(1,2,2)
imagesc(ROI);
ca = gca;
ca.YDir = 'normal';
axis image
shg
