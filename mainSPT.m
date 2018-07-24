% this script is a work in progress, the plan is to fix here the quick
% callibration of the channels. Basically I want to be able to lead the
% movies and built a dataset of [xDim,yDim,channel,frame] where all
% channels represent the same region of the image to a pixel resolution.
clear
close all
clc

% ad general paths that can be usefull
addpath(genpath('Ext'));

% path to the callibration

path2zCal = '..\data\Multiplane\ZCalibration\BeadsZCalibration_1';
path2File = '..\data\Multiplane\TL\TL-OD2-200msExposure_1';
path2Cal = '..\data\Multiplane\PlaneCalib\BeadsCalibrationZStack_1';

detectParam.delta = 6;
detectParam.chi2 = 80;

%% create a Movie Object
mov1 = Core.Movie(path2File);

%% showFrame
mov1.showFrame(22);
%% Calib
calib = Core.calib2D(path2Cal);

calib.calc;

calib.showFrame(22)

%% Calibrate
mpMov = Core.mpMovie(path2File,calib.getCal);

mpMov.calibrate;

mpMov.showFrame(15);

%%
mpLocMov = Core.mpLocMovie(path2File,calib.getCal);
mpLocMov.giveInfo;
mpLocMov.calibrate;
%% getCandidatePos

mpLocMov.findCandidatePos(detectParam);
candidate = mpLocMov.getCandidatePos(24);
mpLocMov.showCandidate(24);

%%
mpLocMov.superResConsolidate;
%% ZCal
zCalMov = Core.zCalMovie(path2zCal,calib.getCal);
%% CandidatePos
zCalMov.giveInfo;
zCalMov.findCandidatePos(detectParam);

zCalMov.superResConsolidate;
%mov1.superResConsolidate(6);

%% showParticles
zCalMov.showParticles(24);

%% ZCalibration
traces =  zCalMov.trackInZ;

%% Show traces
zCalMov.showParticlesTracked(30);%ips
%% ZCalibrate
[zData] = zCalMov.zCalibrate;

%% Show ZCalibration
zCalMov.showZCalibration;

%% Get 3D traces
traces = zCalMov.get3DTraces;

%% Show Traces

zCalMov.showTraces

%% SuperResCalMovie
trackParam.euDistPx = 1; 
trackParam.ellip = 5;
SRCalMov = Core.superResCalMov(path2Cal,calib.getCal);

SRCalMov.giveInfo;
SRCalMov.findCandidatePos(detectParam);

SRCalMov.showCandidate(50);

SRCalMov.superResConsolidate;

SRCalMov.showParticles(50);
%% SuperResCalibrate

SRCalMov.superResCalibrate(trackParam);

%% example of a frame list I will grow this into the frame object
frameList = mcodekit.list.dl_list();
for i = 1:8
    tmp = data(:,:,i,1);
    imP = Core.imPlane(tmp);
    imP.setPixSizeNm(100);
    imP.setTime(uint16(1)); 
    frameList.append_key(imP);   
end

%% detect particles
% for a water immersion obj the best-fit gasuss to the PSF has 
% sigma = 0.25 wavelength / NA
objNA  = 1.2;
emWave = 600;
sigma_nm = 0.25 * emWave/objNA;
FWHMnm = sigma_nm * sqrt(8*log(2));         

GLRTprops.delta  = 6;
GLRTprops.pxSnm  = 100;
GLRTprops.FWHMnm = FWHMnm;
GLRTprops.chi2   = 80;

rTh = 5; % in pixels
ROIrad = 10;

pList = [];
p = [];
partList = mcodekit.list.dl_list();

for fIdx = frame(1:10)
    
    sfData = data(:,:,:,fIdx);
    imSize = size(sfData);
    % detect particles
    [consLoc,totLoc] = mpSetup.localize(sfData, rTh, GLRTprops);
    % build ROIs
    [ROIs] = Misc.getROIs(consLoc,ROIrad,imSize);
    
    for i = 1:size(consLoc,1)
        tmpLoc = consLoc(i,:);
        tmpROI = ROIs(i,:);
        p = Core.particle(tmpLoc,fIdx,tmpROI,sfData);

        pList = [pList, p];
        partList.append_key(p);
        
    end
    disp(['done for frame ' num2str(fIdx)])

end

%%
objNA    = 1.2;
emWave   = 600;
pxSizeNm = 100;
tic
pList.setPSFprops(objNA, emWave, pxSizeNm);
toc

tic
iterator = partList.get_iterator(); 
            
while (iterator.has_next())

    pTmp = iterator.next();
    pTmp.setPSFprops(objNA, emWave, pxSizeNm);
%     idx = idx + 1;

end
toc

%%
ptest = partList.get_key(10);
%%

tic
pList.superResolve;
disp('Done with SR-loc')
toc

tic
iterator = partList.get_iterator(); 
            
while (iterator.has_next())

    pTmp = iterator.next();
    pTmp.superResolve();
%     idx = idx + 1;

end
toc

ptest1 = pList(2);
ptest2 = partList.get_key(2);

%%
tic
test = findobj(pList,'frame',1);
toc

tic
iterator = partList.get_iterator(); 
ftotal = zeros(partList.size_,1);
idx = 0;
while (iterator.has_next())

    idx = idx + 1;
    pTmp = iterator.next();
    ftotal(idx) = pTmp.frame;

end

idxList = find(ftotal==1);
toc
bla = [];
for idx = idxList'
    idx
    tmp = partList.get_key(idx);
    bla = [bla tmp];
end
%%
% test = findobj(pList,'frame',1);
% tmpVal = cat(1,test.superResLoc);
% tmpX = tmpVal(:,1);
% tmpY = tmpVal(:,2);
% tmpZ = tmpVal(:,3);
% % scatter3 (tmpX,tmpY, tmpZ)
% scatter (tmpX,tmpY,'kx')
% 
% % axis image
% shg

cols = {'k','r','g','b','y'};

for ii = 1:250
    
    test = findobj(pList,'frame',ii);
    cidx = ceil(ii/50);
    tmpVal = cat(1,test.superResLoc);
    
    tmpX = tmpVal(:,1);
    tmpY = tmpVal(:,2);
    tmpZ = tmpVal(:,3);
%     scatter3 (tmpX,tmpY, tmpZ)
    subplot(1,2,1)
    scatter (tmpX,tmpY,[cols{cidx} 'x'])
    hold on
    subplot(1,2,2)
    scatter (tmpX,tmpY,[cols{cidx} 'x'])
    hold on
end
subplot(1,2,1)
hold off
axis image
d = .6;
c = [270.5,243.6];
xlim([c(1)-d c(1)+d])
ylim([c(2)-d c(2)+d])
% ylim([242 243])

subplot(1,2,2)
hold off
axis image
d = .6;
c = [430.3,277.2];
xlim([c(1)-d c(1)+d])
ylim([c(2)-d c(2)+d])
shg
%%
% subplot(1,2,1)
% xlim([269.5,271.5])
% ylim([243, 244])

% make a common list of sm detections? should I have a test for seeing a
% molecule in at least 3 (or X) planes? once I have the common list I have
% to cropt the ROIs [dx, dy, 8] and do the fine localization. We will for
% the moment pick the best as I do not have a super-resolved registration
% matrix between all channels. This should be generated in order to use
% info from multiple planes to increase fit accuracy. it is from z tacks of
% beads that we can create such a registration.
