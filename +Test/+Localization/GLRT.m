function GLRT()
doPlot = true;
pxSize = 95; % in nm
imSize = 200; % in px
setupPSFWidth = 220; %in nm (Calculated in Focus using PSFE plate, on the
%15/02/2018 Exc wavelength = 532nm;

prompt = {'Enter number of frame to simulate: ',...
    'Enter number of molecules/frame: ',...
    'Enter a type of noise to add (none, Gaussian or Poisson):',...
    'Enter Signal to noise ratio (for Gaussian): ',...
    'Enter background:','Enter emitter intensity:',...
    'Enter emitter intensity distribution width:',...
    'Enter emitter max ellipticity:'};

dlgTitle = 'Simulation Parameters input';
numLines = 1;
defaultVal = {'500','10','Gaussian','10','500','10000','100','1'};
answer = inputdlg(prompt, dlgTitle,numLines,defaultVal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Check USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(~isempty(answer),'User canceled input dialog, Simulation was aborted')
nSim = str2double(answer(1));
assert(~isnan(nSim),'Number of simulation should be numerical');%If not a number
%expressed in string str2double yield NaN, isnumeric yield true on NaN so
%isnan is the only way to check.
nEm = str2double(answer{2});

noiseType = answer(3);
noiseType = noiseType{1};
assert(isnan(str2double(noiseType)),'Type of noise should not be a number');

S2N = str2double(answer(4));
assert(~isnan(S2N),'Variance should be numerical');

bkgMean = str2double(answer(5));
assert(~isnan(bkgMean),'Background should be numerical');

maxCount = str2double(answer(6));
assert(~isnan(maxCount),'Max count should be numerical');

emIntSigma = str2double(answer(7));
assert(~isnan(emIntSigma),'Max count should be numerical');

emMaxSigma = str2double(answer(8));
assert(~isnan(emMaxSigma),' emitter width should be numerical');
%%%%%%%%%%%%%%%%%%%%%%% END CHECK USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
detector.xSize = imSize;
detector.pxSize = pxSize; %[nm/pix]

emitter.num = nEm;
emitter.meanInt  = maxCount;
emitter.intSigma = emIntSigma;
emitter.FWHM_nm  = setupPSFWidth;
emitter.noiseType= noiseType;
emitter.posRange = [1+round(0.1*imSize) imSize-round(0.1*imSize)];
emitter.maxSigma = emMaxSigma;
bkg.mean = bkgMean;
bkg.SNR  = S2N;

[imStack,simPos,simElip] = EmitterSim.simulateImages(nSim,emitter,detector,bkg);

%%
h = waitbar(0,'Localization...');
totPos = zeros(size(simPos));
avgPos = zeros(1,size(imStack,3));
avgSimPos = avgPos;

for i = 1 : size(imStack,3)
im_in = double(imStack(:,:,i));

delta = 4;
FWHM_pix = setupPSFWidth / pxSize; %[pix]
% for GLRT
chi2 = 24;
[ pos, ~ ] = Localization.smDetection(im_in, delta, FWHM_pix, chi2 );
totPos(1:size(pos,1),:,i) = pos;
avgPos(1,i) = sqrt(mean(pos(:,1))^2 + mean(pos(:,2))^2);
avgSimPos(1,i) = sqrt(mean(simPos(:,1,i))^2 + mean(simPos(:,2,i))^2);

waitbar(i/size(imStack,3),h,sprintf('Localization... - %d / %d percent achieved',round(i/size(imStack,3)*100), 100));
end

if doPlot
figure(1)

imagesc(im_in(:,:,1))
hold on
colormap('hot');
plot(pos(:,1,1),pos(:,2,1),'b+')
hold off
xlabel('Pixel')
ylabel('Pixel')
title('Examplary simulated frame and results of the localization')
axis image
end

disp('------------------------ TEST RESULTS ----------------------------')
fprintf('--------- %d / %d molecules localized \n', size(find(totPos),1),size(find(simPos),1));
fprintf('--------- Center of localization deviated on average from %0.2f pixels\n\n',mean(abs(avgPos-avgSimPos)));

close(h);
end




