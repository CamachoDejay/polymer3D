%Aim of the code:
%Test out the whole procedure that shall be performed in data analysis:
% A) Receiving a stack of images (here simulated)
% B) Localizing the emitters
% C) Determining their position via gradientFit
%We will then check how the results look like for different signal to noise
%ratio and check if gradientFit performs better than the centroid method

%% User Input
nImages  = 1;
% information about emitters
emitter.num      = 10; %Makes at max 10000 fit ==> similar to previous simulation
emitter.meanInt  = 10000;
emitter.intSigma = 100;
emitter.FWHM_nm  = 350;

% Information about detector
detector.xSize  = 100;
detector.pxSize = 105; %[nm/pix]

% information about normal bg
bkg.mean = 1000;
bkg.SNR  = 100;

nSim = nImages*emitter.num;

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