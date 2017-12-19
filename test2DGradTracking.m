%%
%Description: test2DGradTracking, simulate perfect/noisy PSF and use
%gradient fitting to find the center. A table (simResults) is output which
%contain the real value as well as the simulation parameters and the fit 
%value. 

%1) The User will be shown a prompt where he can input some parameter for the
%simulation:
%- Number of PSF to simulate
%- Type of noise to add (none,Gaussian, Poisson).
%- Signal to noise (for Gaussian noise)
%- Backgroud to add
%- Max count (= Gaussian amplitude)

%2) PSF is simulated and noise is added according to the user input

%3) Gradient fit is performed

%4) Output table is generated

%TO DO : Add sigmaX and sigmaY in the table
%TO DO : Compare Sigmas with elipticity from Values 

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doPlot = true;
pix_size = 100; % in nm
im_size = 13; % in px

prompt = {'Enter number of simulation: ',...
    'Enter a type of noise to add (none, Gaussian or Poisson):',...
    'Enter Signal to noise ratio (for Gaussian): ',...
    'Enter background:(1:1000)):','Enter max count:'};
dlgTitle = 'Simulation Parameters input';
numLines = 1;
defaultVal = {'1','none','10','10','100'};
answer = inputdlg(prompt, dlgTitle,numLines,defaultVal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Check USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(~isempty(answer),'User canceled input dialog, Simulation was aborted')
nSim = str2double(answer(1));
assert(~isnan(nSim),'Number of simulation should be numerical');%If not a number
%expressed in string str2double yield NaN, isnumeric yield true on NaN so
%isnan is the only way to check.

noiseType = answer(2);
noiseType = noiseType{1};
assert(isnan(str2double(noiseType)),'Type of noise should not be a number');

S2N = str2double(answer(3));
assert(~isnan(S2N),'Variance should be numerical');

bkg = str2double(answer(4));
assert(~isnan(bkg),'Background should be numerical');
maxCount = str2double(answer(5));
assert(~isnan(maxCount),'Max count should be numerical');

%%%%%%%%%%%%%%%%%%%%%%%%%% END Check USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation

% Allocate memory for storing results
simResults = table(zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),cell(nSim,1), zeros(nSim,1),...
    'VariableNames',{'realX','fitX', 'realY', 'fitY','realElip','cElip','fitElip',...
    'signal2Bkg','signal2Noise','absErrorX','absErrorY','errorFitElip',...
    'noiseType','background'});

noiseProp = struct('S2N',S2N,'bkg',bkg,'maxCount',maxCount);
%simulate data, analyze and store results
for i = 1: nSim
    
    pos_real = [400 + 400*rand(1),400+400*rand(1)];%random number between 1 and 2. (1 yields
    %px No 5 while 2 give pixel number 9 ==> center pixel +-2.
    
    sigX = 200+rand(1)*200;%Generate random number between 0.6 and 1.2
    sigY = 200+rand(1)*200;
    
    xid = 0:im_size-1;
    yid = 0:im_size-1;
    
    xVal = xid.*pix_size;
    yVal = yid.*pix_size;
    pos_pix = (pos_real./pix_size) + 1;
    
    sig = [sigX,sigY];
    ROI = GradientFit.gaus2D(pos_real,sig,xVal,yVal,noiseProp.maxCount); %Generate 2D gaussian
    
    % ROI coor is always the center position
    ROI_coor = [median(1:size(ROI,1)),median(1:size(ROI,1))];
    
    GraR = 4; % The radius of Gradient used for caculation
    
    simResults.signal2Bkg(i) = noiseProp.maxCount/bkg;
    simResults.signal2Noise(i)  = noiseProp.S2N;
    
    % Adding noise onto the "perfect" gaussian
    ROI = GradientFit.generateNoise(ROI,noiseType,noiseProp);
    
    % Do gradient fitting
    [x,y,e,centOut] = Localization.gradFit(ROI,GraR);
   
    xc = (ROI_coor(1) + x);%in px
    yc = (ROI_coor(2) + y);%in px

    %Store the results
    simResults.realX(i)      = (pos_real(1)/pix_size)+1;
    simResults.realY(i)      = (pos_real(2)/pix_size)+1;
    simResults.fitX(i)       = xc;
    simResults.fitY(i)       = yc;
    simResults.realElip(i)   = sigY/sigX;
    simResults.cElip(i)      = centOut.e;
    simResults.fitElip(i)    = e;
    simResults.noiseType(i)  = {noiseType};
    simResults.background(i) = bkg;
    
    
end
%Calculate errors and store them
simResults.absErrorX    = abs(simResults.fitX-simResults.realX);
simResults.absErrorY    = abs(simResults.fitY-simResults.realY);
simResults.errorFitElip = simResults.fitElip-simResults.realElip;

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
