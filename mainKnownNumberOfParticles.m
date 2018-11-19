clear 
clc
close all
%% User input
chi2 = 80;
FWHM_pix = 4;
delta = 20;

%% Create the Movie object

path2File= 'E:\Data\Leuven Data\2018\ZHao\TestCode\400 nm AuNPs 1064 nm laser stepwise - 2 polarization_1';
info.type = 'transmission';
myMov = Core.Movie(path2File,info);
myMov.giveInfo;

%%
myMov.showFrame(1,2)

%% get the data from the movie
fullStack = myMov.getFrame;

%% inversion of the scale cam1 extraction

fullStack = imcomplement(fullStack.Cam1);
%% detection of the center of the beads
firstFrame = double(fullStack(:,:,1)); %Conversion to double for detection

[ pos, ~, ~] = Localization.smDetection(firstFrame, delta, FWHM_pix, chi2 );

pos = round(mean(pos,1));
%% Cropping Movie

fullStack = fullStack(pos(1)-delta:pos(1)+delta, pos(2)-delta:pos(2)+delta,:);
%%
nFrames = size(fullStack,3);
for i = 1:nFrames
    
    %Code for fitting gaussian should be here 
    
end

