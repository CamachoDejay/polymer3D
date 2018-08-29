%The aim of this code is to calculate, from a zStack of a USAF-like target,
%the pixel Size of the setup. It is based on Susana target which was
%showing peaks for the plane best in focus
%a 

clear
close all
clc
%% User input
knownDistance = 10;%in µm
 
%%
fPath = 'E:\Data\Leuven Data\2018\08-Aug\29\PlaneCalibration_25';
fName = 'PlaneCalibration_25_MMStack_Pos0.ome.tif';

fPath = [fPath filesep fName];

% Calculate calibration

[cal] = mpSetup.cali.calculate(fPath, false);

[data] = mpSetup.loadAndCal( fPath, cal );
%% calculate magnification
pxSize = zeros(1,size(data,3));
 figure
 hold on
for i = 1 : size(data,3)
%i = 1   
    idx = cal.neworder(i);
    frameIdx = cal.inFocus(idx).frame;
    
    planeFocus = squeeze(data(:,:,i,frameIdx));
    
    pattern = mean(planeFocus,2);
    % remove background:
    pattern(pattern<3*median(pattern)) = 0;
    
    %find zeros
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                    % Returns Zero-Crossing Indices Of Argument Vector
    zx = zci(pattern);  
    %find where value change from zeros to non zeros
    test = diff(zx);
    [val,peak] = findpeaks(test);
    %corr peak
    corr = cumsum(val-1);
    peak(2:end) = peak(2:end)+corr(1:end-1);
    nMax = length(peak)-1;
    peakLoc = zeros(1,length(peak));
    %search max of the local non zero region
    for j = 1 : nMax+1
        [~,idx] = max(pattern(peak(j):peak(j)+val(j)));
        peakLoc(j) = peak(j)+idx-1;
    end
    
    plot(pattern)

    distanceInPx = max(peakLoc)-min(peakLoc);
    distanceInUm = nMax*knownDistance;
   
    pxSize(i) = distanceInUm*1000/distanceInPx;
    
end
%% plotting
figure
plot(pxSize)

