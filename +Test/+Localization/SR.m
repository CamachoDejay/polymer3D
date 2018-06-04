function SR()
%Description: SRLocalization, simulate perfect/noisy PSF and use the chosen
%algorithm to fit it. A table (simResults) is output which
%contain the real value as well as the simulation parameters and the fit
%value.

%1)INPUT : The User will be shown a prompt where he can input some parameter for the
%simulation:
%- Number of PSF to simulate
%- Type of noise to add (none,Gaussian, Poisson).
%- Signal to noise (for Gaussian noise)
%- Backgroud to add
%- Max count (= Gaussian amplitude)
%- Fitting method to use (Phasor or Gradient)

%2) PSF is simulated and noise is added according to the user input

%3) Fit is performed

%4) OUTPUT: table is generated containing simulated data and fitted results
%as well as user input parameters
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doPlot = true;
pxSize = 95; % in nm
imSize = 13; % in px
setupPSFWidth = 300; %in nm (Calculated in Focus using PSFE plate, on the
%15/02/2018 Exc wavelength = 532nm;

prompt = {'Enter number of simulations: ',...
    'Enter a type of noise to add (none, Gaussian or Poisson):',...
    'Enter Signal to noise ratio (for Gaussian): ',...
    'Enter background:','Enter max count:','Enter the maximum Ellipticity:',...
    'Enter Min pos', 'Enter Max pos',...
    'Enter the method to be used (Phasor or Gradient):'};
dlgTitle = 'Simulation Parameters input';
numLines = 1;
defaultVal = {'10000','Gaussian','8','1000','10000','3','5','9','Phasor'};
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

bckg = str2double(answer(4));
assert(~isnan(bckg),'Background should be numerical');

maxCount = str2double(answer(5));
assert(~isnan(maxCount),'Max count should be numerical');

emMaxSigma = str2double(answer(6));
assert(~isnan(emMaxSigma),'Emitter width should be numerical');

minPos = str2double(answer(7));
assert(~isnan(minPos),'Min position should be numerical');

maxPos = str2double(answer(8));
assert(~isnan(maxPos),'Max position should be numerical');

assert(minPos<= maxPos,'Min position should be smaller than max position');

fitting = answer(9);
assert(or(or(strcmp(fitting,'Phasor'),strcmp(fitting,'phasor')),...
    or(strcmp(fitting,'Gradient'),strcmp(fitting,'gradient'))),...
    'The requested algorithm is unknown, please check spelling');
%%%%%%%%%%%%%%%%%%%%%%%%%% END Check USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation

% Allocate memory for storing results
simResults = table(zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    cell(nSim,1), zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),cell(nSim,1),...
    'VariableNames',{'realX','fitX', 'realY', 'fitY','realElip','cElip','fitElip',...
    'signal2Bkg','signal2Noise','fAbsErrorX','fAbsErrorY','cAbsErrorX',...
    'cAbsErrorY','errorFitElip','noiseType','background','cFitX',...
    'cFitY','fittingMethod'});

noiseProp = struct('S2N',S2N,'bkg',bckg,'maxCount',maxCount);

GraR = 4; % The radius of Gradient used for caculation

%Info for simulation
emitter.num = 1;
emitter.meanInt  = maxCount;
emitter.intSigma = 0.1 * maxCount;
emitter.FWHM_nm  = setupPSFWidth;
emitter.posRange = [minPos maxPos];
emitter.noiseType= noiseType;
emitter.maxSigma = emMaxSigma;
detector.xSize   = imSize;
detector.pxSize  = pxSize; %[nm/pix]
bkg.mean = bckg;
bkg.SNR = S2N;

h = waitbar(0, 'Simulations and Fitting...');
for i = 1: nSim    
    [ROI,simPos,simElip] = EmitterSim.simulateImages(1,emitter,detector,bkg);
    
    % ROI coor is always the center position
    ROI_coor = [median(1:size(ROI,1)),median(1:size(ROI,1))];

    simResults.signal2Bkg(i) = noiseProp.maxCount/bckg;
    simResults.signal2Noise(i)  = noiseProp.S2N; 
 
    if or(strcmp(fitting,'gradient'),strcmp(fitting,'Gradient'))
      [x,y,e,centOut] = Localization.gradFit(ROI,GraR);

    else
     [x,y,e] = Localization.phasor(ROI);
     centOut.x = NaN;
     centOut.y = NaN;
     centOut.e = NaN;
    end
     
    %Test fitting output
    if abs(x) > GraR || abs(y) > GraR || e<=0
        x = NaN;
        y = NaN;
        e = NaN;
    end
    
    xc = (ROI_coor(1) + x);%in px
    yc = (ROI_coor(2) + y);%in px
    
    %Store the results
    simResults.realX(i)      = simPos(1);
    simResults.realY(i)      = simPos(2);
    simResults.fitX(i)       = xc;
    simResults.fitY(i)       = yc;
    simResults.cFitX(i)      = (ROI_coor(1) + centOut(1).x);
    simResults.cFitY(i)      = (ROI_coor(2) + centOut(1).y);
    simResults.realElip(i)   = simElip;
    simResults.cElip(i)      = centOut.e;
    simResults.fitElip(i)    = e;
    simResults.noiseType(i)  = {noiseType};
    simResults.fittingMethod(i) = {fitting};
    simResults.background(i) = bckg;
    simResults.cAbsErrorX(i) = (ROI_coor(1)+centOut(1).x)-simResults.realX(i);
    simResults.cAbsErrorY(i) = (ROI_coor(2)+centOut(1).y)-simResults.realY(i);
    
    waitbar(i/nSim,h, sprintf('Simulation and fitting... %d /100 percent done',round(100*i/nSim)));
end
%Calculate errors and store them
simResults.fAbsErrorX    = simResults.fitX-simResults.realX;
simResults.fAbsErrorY    = simResults.fitY-simResults.realY;
simResults.errorFitElip = simResults.fitElip-simResults.realElip;

%% Plotting and displaying test results

disp('------------------------ TEST RESULTS ----------------------------')
fprintf('--------- %d / %d molecules succesfully localized \n', size(find(~isnan(simResults.fitX)),1),nSim);
fprintf('--------- Center of localization deviated on average from %0.4f pixels in X\n',nanmedian(abs(simResults.fAbsErrorX)));
fprintf('--------- Center of localization deviated on average from %0.4f pixels in Y\n',nanmedian(abs(simResults.fAbsErrorY)));
fprintf('--------- Ellipticity deviated on average from %0.2f\n\n',nanmedian(abs(simResults.errorFitElip)));

if doPlot
    figure(1)
    subplot(1,3,1)
    imagesc(ROI);
    xlabel('x-pos')
    ylabel('y-pos')
    axis image
    title('Example of simulated ROI')
    
    subplot(1,3,2)
    histogram(simResults.fAbsErrorX)
    hold on
    histogram(simResults.fAbsErrorY)
    xlabel('Error (pixel)')
    ylabel('Occurrence')
    
    subplot(1,3,3)
    scatter(simResults.realElip,simResults.fitElip)
    xlabel('Simulated ellipticity')
    ylabel('Fitted ellipticity')
end

close(h)

end