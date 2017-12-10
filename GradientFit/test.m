clear 
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noise = 'Gaussian'; %'Gaussian', 'Poisson','both?'
doPlot = false;
nSim = 1;
pix_size = 0.25;
im_size = 13; % in px
noiseFactor = 10;%For Gaussian

%%%%%%%%%%%%%%%%%%%%%%%%%%%% END USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocate memory for storing results
simResults = table(zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    zeros(nSim,1), zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    'VariableNames',{'realX','fitX', 'realY', 'fitY','elipticity','absErrorX',...
    'absErrorY','relErrorX','relErrorY'}); 

for i = 1: nSim

    pos_real = [1 + rand(1),1+rand(1)];%random number between 1 and 2. (1 yields
    %px No 5 while 2 give pixel number 9 ==> center pixel +-2.

    sigX = 0.6+rand(1)*0.6;%Generate random number between 0.6 and 1.2
    sigY = 0.6+rand(1)*0.6;

    xid = 0:im_size-1;
    yid = 0:im_size-1;

    xVal = xid.*pix_size;
    yVal = yid.*pix_size;
    pos_pix = (pos_real./pix_size) + 1;

    sig = [sigX,sigY];
    ROI = gaus2D(pos_real,sig,xVal,yVal); %Generate 2D gaussian
    
    %Casting to integers
   
    
    % ROI coor is always the center position
    ROI_coor = [median(1:size(ROI,1)),median(1:size(ROI,1))];

    GraR = 4; % The radius of Gradient used for caculation

    % Adding noise onto the "perfect" gaussian
    switch noise
        case 'Gaussian'
            bkg = 100;
            ROI = ROI +bkg;
            ROI = ROI + 10* randn(size(ROI));
           % ROI = uint16(ROI);
            %add Gaussian distributed noise
        case 'Poisson'
            %add Poisson distributed noise
        case 'both'
            %add Gaussian and poisson noise
        otherwise
           %  ROI = uint16(ROI);
    end

    [x,y,e] = GradientFit(ROI,GraR);% Do gradient fitting
    
    xc = (ROI_coor(1) + x);
    yc = (ROI_coor(2) + y);
    
    simResults.realX(i) = pos_real(1);
    simResults.realY(i) = pos_real(2);
    simResults.fitX(i)  = (xc-1)*pix_size;
    simResults.fitY(i)  = (yc-1)*pix_size;
    simResults.elipticity(i) = e;
   
    

end
 simResults.absErrorX = abs(simResults.fitX-simResults.realX);
 simResults.absErrorY = abs(simResults.fitY-simResults.realY);
 simResults.relErrorX = simResults.absErrorX./simResults.realX;
 simResults.relErrorY = simResults.absErrorY./simResults.realY;
 
%to Display in the command line
%fprintf('x pos [pix]: %.4g  [fit]: %.4g [real]: %.4g \t \n',xc,(xc-1)*pix_size,pos_real(1))
%fprintf('y pos [pix]: %.4g  [fit]: %.4g [real]: %.4g \t \n',yc,(yc-1)*pix_size,pos_real(2))
%fprintf('elipt: %.4g \n',e)
% fprintf('x out: %.4g \n',x)
% fprintf('y out: %.4g \n',y)


if doPlot
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
end
