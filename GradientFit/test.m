%%
clear 
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doPlot = true;
pix_size = 0.25;
im_size = 13; % in px


prompt = {'Enter number of simulation: ',...
    'Enter a type of noise to add (none, Gaussian or Poisson):',...
    'Enter variance to use for Gaussian noise: ',...
    'Enter background:(1:1000)):'};
dlgTitle = 'Simulation Parameters input';
numLines = 1;
defaultVal = {'1','none','0.1','10'};
answer = inputdlg(prompt, dlgTitle,numLines,defaultVal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%% Check USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(answer)
    warning('User canceled input dialog, Simulation will run with default input');
    nSim  = 1;
    noise = 'none';
    gVar  = 0.1;
    bkg   = 10;
else

nSim = str2double(answer(1));
assert(~isnan(nSim),'Number of simulation should be numerical');%If not a number
%expressed in string str2double yield NaN, isnumeric yield true on NaN so
%isnan is the only way to check.

noise = answer(2);
noise = noise{1};
assert(isnan(str2double(noise)),'Type of noise should not be a number'); 

gVar = str2double(answer(3));
assert(~isnan(gVar),'Variance should be numerical');

bkg = str2double(answer(4));
assert(~isnan(bkg),'Background should be numerical');
end

%%%%%%%%%%%%%%%%%%%%%%%%%% END Check USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation

% Allocate memory for storing results
simResults = table(zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    zeros(nSim,1), zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    cell(nSim,1), zeros(nSim,1), zeros(nSim,1),...
    'VariableNames',{'realX','fitX', 'realY', 'fitY','elipticity','absErrorX',...
    'absErrorY','relErrorX','relErrorY','noiseType','background','gVariance'}); 

%simulate data, analyze and store results 
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
   
    
    % ROI coor is always the center position
    ROI_coor = [median(1:size(ROI,1)),median(1:size(ROI,1))];

    GraR = 4; % The radius of Gradient used for caculation

    % Adding noise onto the "perfect" gaussian
    switch noise
        case 'Gaussian'
            %Add background to avoid negative value after adding noise           
            tmpROI = ROI +bkg;%to avoid negative value with Gaussian noise
            
            gNoise = sqrt(gVar).*randn(size(tmpROI))+ 0;%randn yield var = 1 multiplying
            %the noise will multiply variance the same amount to the power
            %of 2, thus we sqrt the variance demanded by the user.
            ROI = round (tmpROI+gNoise);%Rounded to obtain integers while keeping
            %double type for gradient function.
            simResults.gVariance(i) = gVar;
        case 'Poisson'
           tmpROI = uint16(ROI +bkg); %Cast to int for imnoise
           ROI = double(round(imnoise(tmpROI,'poisson')));   
           simResults.gVariance = 0;
      %   case 'both'
          %ROI = ROI +bkg;
          % ROI = imnoise(ROI,'Gaussian',0,gVar);
          % ROI = imnoise(ROI,'poisson');
          
        otherwise
           ROI = round(ROI);%Rounded to obtain integers while keeping
            %double type for gradient function.
            
           simResults.gVariance = 0;
    end

    [x,y,e] = GradientFit(ROI,GraR);% Do gradient fitting
    
    xc = (ROI_coor(1) + x);%in px
    yc = (ROI_coor(2) + y);%in px
    
    %Store the results 
    simResults.realX(i)      = pos_real(1);
    simResults.realY(i)      = pos_real(2);
    simResults.fitX(i)       = (xc-1)*pix_size;
    simResults.fitY(i)       = (yc-1)*pix_size;
    simResults.elipticity(i) = e;
    simResults.noiseType(i)  = {noise};
    simResults.background(i) = bkg; 
   
end
%Calculate errors and store them
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
