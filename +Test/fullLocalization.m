function fullLocalization()
%Aim of the code:
%Test out the whole procedure that shall be performed in data analysis:
% A) Receiving a stack of images (here simulated)
% B) Localizing the emitters
% C) Determining their position via gradientFit
%We will then check how the results look like for different signal to noise
%ratio and check if gradientFit performs better than the centroid method

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doPlot = true;
pxSize = 95; % in nm
imSize = 13; % in px
setupPSFWidth = 300; %in nm (Calculated in Focus using PSFE plate, on the
%15/02/2018 Exc wavelength = 532nm;

prompt = {'Enter number of simulation: ',...
    'Enter a type of noise to add (none, Gaussian or Poisson):',...
    'Enter Signal to noise ratio (for Gaussian): ',...
    'Enter background:','Enter max count:', 'Enter Min pos',...
    'Enter Max pos', 'Enter the method to be used (Phasor or Gradient):'};
dlgTitle = 'Simulation Parameters input';
numLines = 1;
defaultVal = {'10000','Gaussian','8','500','10000','5','9','Phasor'};
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

minPos = str2double(answer(6));
assert(~isnan(minPos),'Min position should be numerical');

maxPos = str2double(answer(7));
assert(~isnan(maxPos),'Max position should be numerical');

assert(minPos<= maxPos,'Min position should be smaller than max position');

fitting = answer(8);
assert(or(or(strcmp(fitting,'Phasor'),strcmp(fitting,'phasor')),...
    or(strcmp(fitting,'Gradient'),strcmp(fitting,'gradient'))),...
    'The requested algorithm is unknown, please check spelling');
%%%%%%%%%%%%%%%%%%%%%%%%%% END Check USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Allocation memory for results storage
% Allocate memory for storing results
simResults = table(zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),zeros(nSim,1),...
    'VariableNames',{'realX','fitX', 'realY', 'fitY','realElip','cElip',...
    'fitElip','signal2Noise','fAbsErrorX','fAbsErrorY','cAbsErrorX',...
    'cAbsErrorY','errorFitElip'});

%% Image stack simulation 

%Let us generate a stack of images
[imStack,simPos,simElip] = Misc.simulateImages(nImages,emitter,detector,bkg);

simResults.realX(:) = simPos(:,1);
simResults.realY(:) = simPos(:,2);
simResults.realElip(:) = simElip(:,1);
simResults.signal2Noise(:) = bkg.SNR;


%% Localization and gradient Fitting

szWindow = 6;
imBg = zeros(2*szWindow+1,2*szWindow+1);
imBg = imBg+bkg.mean;
xSize = size(imStack,2);
ySize = size(imStack,1);
GraR = 2; % The radius of Gradient used for caculation
countLocMol = 0;
tic
for i=1:size(imStack,3)
    
    im_in = double(imStack(:,:,i));
    delta = 4;
    FWHM_pix = emitter.FWHM_nm / detector.pxSize; %[pix]
    % for GLRT
    chi2 = 24;
    %Localzation occurs here
    [ pos, inten ] = Localization.smDetection(im_in, delta, FWHM_pix, chi2 );
    countLocMol = countLocMol + size(pos,1);
    for j=1:size(pos,1)
        %Extract a roi around the localized emitter
        [ roi_lims ] = EmitterSim.getROI(pos(j,1), pos(j,2), szWindow, xSize, ySize);
        % im(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2)) = ...
          %            im(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2)) + im_roi;
       
        
        %Check if full Roi could be taken
        if(abs(roi_lims(1)-roi_lims(2))~=2*szWindow || abs(roi_lims(3)-roi_lims(4))~=2*szWindow)
            %if not give some padding
%             padX = 2*szWindow - (roi_lims(2)-roi_lims(1));
%             padY = 2*szWindow - (roi_lims(4)-roi_lims(3));
            ROI = imBg;
            ROI(1:abs(roi_lims(3)-roi_lims(4))+1,1:abs(roi_lims(1)-roi_lims(2))+1) =...
                imStack(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2),i);
            
        else
            ROI = imStack(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2),i);
%             padX = 0;
%             padY = 0;
        end 
              
        % Do gradient fitting
        [x,y,e,centOut] = Localization.gradFit(ROI,GraR);
        
        % Position in the initial referential
        xc = (round(pos(j,1)) + x);%in px
        yc = (round(pos(j,2)) + y);
        centX = (round(pos(j,1)) +centOut.x);
        centY = (round(pos(j,2)) +centOut.y);
        
        %Check to which initially simulated molecule the value correspond
        [row,~] = find(simResults.realX((i-1)*emitter.num+1:i*emitter.num)<xc+2 & ...
            simResults.realX((i-1)*emitter.num+1:i*emitter.num)>xc-2 & ...
            simResults.realY((i-1)*emitter.num+1:i*emitter.num)<yc+2 & ...
            simResults.realY((i-1)*emitter.num+1:i*emitter.num)>yc-2);
        
        if (length(row)==1)
            simResults.fitX((i-1)*emitter.num+row)    = xc;
            simResults.fitY((i-1)*emitter.num+row)    = yc;
            simResults.fitElip((i-1)*emitter.num+row) = e;
            simResults.cElip((i-1)*emitter.num+row)   = centOut.e;
            simResults.fAbsErrorX((i-1)*emitter.num+row) =...
                xc-simResults.realX((i-1)*emitter.num+row);
            simResults.fAbsErrorY((i-1)*emitter.num+row) =...
                yc-simResults.realY((i-1)*emitter.num+row);
            
            simResults.cAbsErrorX((i-1)*emitter.num+row) =...
                centX-simResults.realX((i-1)*emitter.num+row);
            simResults.cAbsErrorY((i-1)*emitter.num+row) =...
                centY-simResults.realY((i-1)*emitter.num+row);
            
            simResults.errorFitElip((i-1)*emitter.num+row) =...
                e-simResults.realElip((i-1)*emitter.num+row);
        else
             simResults.fitX((i-1)*emitter.num+row)    = NaN;
        end
          
    end

end
toc
end