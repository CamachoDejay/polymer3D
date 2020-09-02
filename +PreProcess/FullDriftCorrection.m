function [corrData] = FullDriftCorrection(data)

%%%%%%%%%%%%%%% DRIFT CORRECTION USING IMAGE CORRELATION %%%%%%%%%%%%%%%%%%
%By Boris Louis, Version of the 30th April 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%USER INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

correlationInfo.corrSz = 100; %in px. Radius of the ROI used for correlation
%correlation function
correlationInfo.driftPeriod = 1; %in Frame, Number of frame that are averaged
%for driftCalculation ==> 1 mean that drift is calculated for each frame
scalingFactor = 1;%Used for interpolation in sub-pixel Drift correction 
%objects of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END USER INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EXPLANATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code aim to apply a drift correction to a bunch of images contained
%in a folder. This is a stand alone code as all the function needed to use
%it are encoded as local function in the very same file. The step of the
%Drift correction are described below:

% #1 User can choose a folder which is then loaded to matlab
% #2 For each SPE file in the folder the following is performed:
%    #2.1 A rough localization of the molecule/crystals is performed
%    #2.2 The image is cropped around the previously found region (limit
%    the computing time by reducing the size of the images to correlate
%    #2.3 The cropped movies are correlated frame by frame to the first
%    frame and the shift between each frame is calculated
%    #2.4 The frame are then shifted with respect to the first frame in
%    order to correct the drift on the movies
%    #2.5 The movies are saved either as .mat files or back to .SPE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END EXPLANATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Drift Correction
%Drift Correction occurs here
[corrData,Drift] = PreProcess.CorrelationDrift(data,...
    scalingFactor,correlationInfo);

end

