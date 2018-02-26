%Description: test2DGradTracking, simulate perfect/noisy PSF and use
%gradient fitting to find the center. A table (simResults) is output which
%contain the real value as well as the simulation parameters and the fit
%value.

%1)INPUT : The User will be shown a prompt where he can input some parameter for the
%simulation:
%- Number of PSF to simulate
%- Type of noise to add (none,Gaussian, Poisson).
%- Signal to noise (for Gaussian noise)
%- Backgroud to add
%- Max count (= Gaussian amplitude)

%2) PSF is simulated and noise is added according to the user input

%3) Gradient fit is performed

%4) OUTPUT: table is generated containing simulated data and fitted results
%as well as user input parameters

clear
close all
clc
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doPlot = true;
pxSize = 100; % in nm
imSize = 21; % in px
filtering = true;
setupPSFWidth = 220; %in nm (Calculated in Focus using PSFE plate, on the
%15/02/2018 Exc wavelength = 532nm;

prompt = {'Enter number of simulation: ',...
    'Enter a type of noise to add (none, Gaussian or Poisson):',...
    'Enter Signal to noise ratio (for Gaussian): ',...
    'Enter background:','Enter max count:', 'Enter Min pos', 'Enter Max pos'};
dlgTitle = 'Simulation Parameters input';
numLines = 1;
defaultVal = {'10000','Gaussian','8','1000','100','5','9'};
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

minPos = str2double(answer(6));
assert(~isnan(minPos),'Min position should be numerical');

maxPos = str2double(answer(7));
assert(~isnan(maxPos),'Max position should be numerical');

assert(minPos<= maxPos,'Min position should be smaller than max position');

%%%%%%%%%%%%%%%%%%%%%%%%%% END Check USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation

% Allocate memory for storing results
simResults = table(zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    cell(nSim,1), zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    'VariableNames',{'realX','fitX', 'realY', 'fitY','realElip','cElip','fitElip',...
    'signal2Bkg','signal2Noise','fAbsErrorX','fAbsErrorY','cAbsErrorX',...
    'cAbsErrorY','errorFitElip','noiseType','background','cFitX',...
    'cFitY'});

noiseProp = struct('S2N',S2N,'bkg',bkg,'maxCount',maxCount);
%simulate data, analyze and store results
halfWay = round(nSim/2);
xid = 0:imSize-1;
yid = 0:imSize-1;

xVal = xid.*pxSize;
yVal = yid.*pxSize;
GraR = 4; % The radius of Gradient used for caculation
tic
for i = 1: nSim
    
    pos_real = [minPos*pxSize-pxSize+ rand(1)*((maxPos-minPos)*pxSize),...
        minPos*pxSize-pxSize+ rand(1)*((maxPos-minPos)*pxSize)];%random number between 400 and 800. (400 yields)
    %     %px No 5 while 800 give pixel number 9 ==> center pixel +-2.
    if i <halfWay
        % pos_real = [600,600];
        sigX = setupPSFWidth+rand(1)*setupPSFWidth*2;
        sigY = setupPSFWidth+rand(1)*10;
        %           sigX = setupPSFWidth;
        %           sigY = setupPSFWidth;
    else
        sigX = setupPSFWidth+rand(1)*10;
        sigY = setupPSFWidth+rand(1)*setupPSFWidth*2;
        %           sigX = setupPSFWidth;
        %           sigY = setupPSFWidth;
    end
    
    pos_pix = (pos_real./pxSize) + 1;
    
    sig = [sigX,sigY];
    ROI = Misc.gaus2D(pos_real,sig,xVal,yVal,noiseProp.maxCount); %Generate 2D gaussian
    
    % ROI coor is always the center position
    ROI_coor = [median(1:size(ROI,1)),median(1:size(ROI,1))];

    simResults.signal2Bkg(i) = noiseProp.maxCount/bkg;
    simResults.signal2Noise(i)  = noiseProp.S2N;
    
    % Adding noise onto the "perfect" gaussian
    ROI = Misc.generateNoise(ROI,noiseType,noiseProp);
    
    if filtering
        %ROI = imgaussfilt(ROI,2);
        ROI = medfilt2(ROI,[2 2],'symmetric');
    end
    % Do gradient fitting
    [x,y,e,centOut] = Localization.gradFit(ROI,GraR);
    
    %Test fitting output
    if abs(x) > GraR || abs(y) > GraR
        x = NaN;
        y = NaN;
        e = NaN;
    end
    
    xc = (ROI_coor(1) + x);%in px
    yc = (ROI_coor(2) + y);%in px
    
    %Store the results
    simResults.realX(i)      = (pos_real(1)/pxSize)+1;
    simResults.realY(i)      = (pos_real(2)/pxSize)+1;
    simResults.fitX(i)       = xc;
    simResults.fitY(i)       = yc;
    simResults.cFitX(i)      = (ROI_coor(1) + centOut(1).x);
    simResults.cFitY(i)      = (ROI_coor(2) + centOut(1).y);
    simResults.realElip(i)   = sigY/sigX;
    simResults.cElip(i)      = centOut.e;
    simResults.fitElip(i)    = e;
    simResults.noiseType(i)  = {noiseType};
    simResults.background(i) = bkg;
    simResults.cAbsErrorX(i) = (ROI_coor(1)+centOut(1).x)-simResults.realX(i);
    simResults.cAbsErrorY(i) = (ROI_coor(2)+centOut(1).y)-simResults.realY(i);
    
end
%Calculate errors and store them
simResults.fAbsErrorX    = simResults.fitX-simResults.realX;
simResults.fAbsErrorY    = simResults.fitY-simResults.realY;
simResults.errorFitElip = simResults.fitElip-simResults.realElip;
toc
disp('Simulation Done');

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
