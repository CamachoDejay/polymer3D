%%%%%%%%%%%%%%%%%%%%%%%%% MAIN TRACKING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
% The aim this code is to plot/view tracking data in a nice way. Coding it
% outside of the object allows for a bit of flexibility.
clear;
close all;
clc;
%% USER INPUT

filePath = 'F:\Data\Leuven Data\2019\04 - April\3\XYZ - CS\X';
dim = 'x';
period = 20;
idx2Plot = 4;
stepApplied = 200;
%% LOADING

fileName = [filePath filesep 'trackResults.mat'];

trackData = load(fileName);
name = fieldnames(trackData);
trackData = trackData.(name{1});

%% Process Data

traces = trackData.traces(:,1);
[stepMot,mot] = getMotor(trackData,dim);
fileRef = cell2mat(trackData.traces(:,2));
nFiles = max(fileRef);
nTraces = length(traces);
lenTraces = cellfun(@height,traces);

[maxLen,idx] = max(lenTraces);
maxFrame = max(traces{idx}.frame);
idx2Mot = [0 period:period:maxFrame+1];

precPerTrace = zeros(nTraces,1);
accPerTrace = zeros(nTraces,1);

precMPerTrace = zeros(nTraces,1);
accMPerTrace = zeros(nTraces,1);
intPerTrace = zeros(nTraces,1);
SNRPerTrace = zeros(nTraces,1);

counter = 0;
precMPerStep = [];
accMPerStep  = [];
precPerStep  = [];
accPerStep   = [];

figure
hold on
allData = nan(nTraces,maxLen);
allDataM = allData;
for i =1:nTraces
    
    currTrace = traces{i};
    currMot   = mot{i}(2:end);
    currStep  = stepMot{i};
    stepSize = max(currStep)*1000;
    lenCTrace = height(currTrace);
    
    data2Plot = getData2Plot(currTrace,dim);
    prec = zeros(length(idx2Mot)-1,1);
    acc = prec;
    precM = zeros(length(idx2Mot)-1,1);
    accM = precM;
    int  = prec;
    SNR  = prec;
    frames = currTrace.frame;

    for j = 1: length(idx2Mot)-1

       idx = idx2Mot(j)+1:idx2Mot(j+1);
       [~,idx2Frame] = intersect(frames,idx);
       if ~isempty(idx2Frame)
           allData(i,idx2Frame) = data2Plot(idx2Frame,1);
           allDataM(i,idx2Frame) = data2Plot(idx2Frame,2);
           prec(j)  = std(data2Plot(idx2Frame,1));
           acc(j)   = mean(data2Plot(idx2Frame,1));
           precM(j) = std(data2Plot(idx2Frame,2));
           accM(j)  = mean(data2Plot(idx2Frame,2));
           int(j)   = mean(currTrace.intensity(idx2Frame));
           SNR(j)   = mean(currTrace.SNR(idx2Frame));

       end
    end
    plot(data2Plot(:,1))
    %clean data
    prec(prec==0)   = [];
    acc(acc==0)     = [];
    precM(precM==0) = [];
    accM(accM==0)   = [];
    int(int==0)     = [];
    SNR(SNR==0)     = [];

    precPerStep  = [precPerStep; prec];
    precMPerStep  = [precMPerStep; precM];
    precPerTrace(i)  = mean(prec);
    precMPerTrace(i)  = mean(precM);
    intPerTrace(i) = mean(int);
    SNRPerTrace(i) = mean(SNR);

    if length(prec)>1

        accPerStep = [accPerStep; abs(diff(acc))-stepSize];
        accMPerStep = [accMPerStep; abs(diff(accM))-stepSize];
        accPerTrace(i) = mean(abs(diff(acc)));
        accMPerTrace(i) = mean(abs(diff(accM)));

    end

    counter = counter+1;
    
end
accPerTrace(accPerTrace==0) = [];
precPerTrace(precPerTrace==0)   = [];
accMPerTrace(accMPerTrace==0) = [];
precMPerTrace(precMPerTrace==0)   = [];

accuracy = nanmean(abs(accPerTrace-stepSize));
precision  = nanmean(precPerTrace);

accuracyM = nanmean(abs(accMPerTrace-stepSize));
precisionM  = nanmean(precMPerTrace);

intensity  = nanmean(intPerTrace);
SNR        = nanmean(SNRPerTrace);


fprintf('The average tracking precision for best focus is %d nm based on %d traces\n',round(precision),counter);
fprintf('The average tracking accuracy for best focus is %d nm  based on %d traces\n',round(accuracy),counter);
fprintf('The FWHM for best focus is %d nm based on %d traces\n',round(2*sqrt(2*log(2))*precision),counter)

fprintf('The average tracking precision for mean is %d nm based on %d traces\n',round(precisionM),counter);
fprintf('The average tracking accuracy for mean is %d nm  based on %d traces\n',round(accuracyM),counter);
fprintf('The FWHM for mean is %d nm based on %d traces\n',round(2*sqrt(2*log(2))*precisionM),counter)

fprintf('The average intensity is %d photons\n',intensity);
fprintf('The average SNR is %d \n',SNR);

figure
subplot(1,2,1)
histogram(precPerStep)
title('Precision for best focus');
subplot(1,2,2)
histogram(accPerStep)
title('Accuracy for best focus')

figure
subplot(1,2,1)
histogram(precMPerStep)
title('Precision for Mean');
subplot(1,2,2)
histogram(accMPerStep)
title('Accuracy for Mean')



%% Plotting
figure 
hold on
%select the traces to plot
idx = fileRef==idx2Plot;
tracePlot = traces(idx);
lenTraces = lenTraces(idx);
idx2 =lenTraces==max(lenTraces);
tracePlot = tracePlot(idx2);

for i = 1:length(tracePlot)
    currTrace = tracePlot{i};
    currMot = mot{i};
    data2Plot = getData2Plot(currTrace,dim);
    scatter(1:length(data2Plot),data2Plot(:,1),'filled')
    
end
motPlot = mot{i}*1000;
if ~strcmp(dim,'z')
    motPlot = abs(motPlot);
end
motPlot = motPlot(currTrace.frame);
motPlot = motPlot - mean(motPlot);
plot(motPlot,'-r','LineWidth',2);

%plot with errorBar
allStd = zeros(nFiles,size(allData,2));
for i = 1:nFiles
    idx = fileRef ==i;
    allStd(i,:) = nanstd(allData(idx,:),1);     
    
end
error = mean(allStd,1);

figure
h= Plotting.shadedErrorBar(1:length(error),motPlot,error);

%% function

function [data2Plot] = getData2Plot(currTrace,dim)

    switch(dim)
        
        case 'x'
            data2Plot(:,1) = currTrace.col-mean(currTrace.col);
            data2Plot(:,2) = currTrace.colM-mean(currTrace.colM);
            
        case 'y'
            data2Plot(:,1) = currTrace.row-mean(currTrace.row);
            data2Plot(:,2) = currTrace.rowM-mean(currTrace.rowM);
        
        case 'z'
            
            data2Plot(:,1) = currTrace.z-mean(currTrace.z);
            data2Plot(:,2) = currTrace.zM-mean(currTrace.zM);
    end

end

function [step,mot] = getMotor(trackData,dim)
    
    switch(dim)
        case 'x'
           
            mot  = trackData.traces(:,4);
            step = trackData.traces(:,3);
        
        case 'y'
        
            mot  = trackData.traces(:,6);
            step = trackData.traces(:,5);
        
        case 'z'
        
            mot  = trackData.traces(:,8) ;
            step = trackData.traces(:,7);
    
    end

end