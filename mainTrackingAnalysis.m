%%%%%%%%%%%%%%%%%%%%%%%%% MAIN TRACKING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
% The aim this code is to plot/view tracking data in a nice way. Coding it
% outside of the object allows for a bit of flexibility.
clear;
close all;
clc;
%% USER INPUT

filePath = 'D:\TmpData\Extended\OD25\Z';
dim = 'z';
period = 20;
idx2Plot = 1;
stepApplied = 200;
%% LOADING

fileName = [filePath filesep 'trackResults.mat'];

trackData = load(fileName);
name = fieldnames(trackData);
trackData = trackData.(name{1});

%% Process Data

traces = trackData.traces(:,1);
[stepMot,mot] = getMotor(trackData,dim);

lenTraces = cellfun(@height,traces);
%only keep traces that have been measured for one period at least
traces(lenTraces<period)= [];
trackData.traces(lenTraces<period,:) = [];
nTraces = length(traces);

fileRef = cell2mat(trackData.traces(:,2));
nFiles = max(fileRef);

[maxLen,idx] = max(lenTraces);
maxFrame = max(traces{idx}.t);
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
allData = nan(nTraces,length(mot{1}));
allDataM = allData;
for i =1:nTraces
    
    currTrace = traces{i};
    currMot   = mot{i}(2:end);
    currStep  = stepMot{i};
    stepSize = max(currStep)*1000;
    lenCTrace = height(currTrace);
    
    data2Plot = getData2Plot(currTrace,dim,mot{i});
    prec = zeros(length(idx2Mot)-1,1);
    acc = prec;
    precM = zeros(length(idx2Mot)-1,1);
    accM = precM;
    int  = prec;
    SNR  = prec;
    frames = currTrace.t;
    %slicing the data in each "step motion"
    for j = 1: length(idx2Mot)-1

       idx = idx2Mot(j)+1:idx2Mot(j+1);
       [~,idx2Frame] = intersect(frames,idx);
       if ~isempty(idx2Frame)
           %to calculate inter particle std we remove the mean
           allData(i,frames(idx2Frame)) = data2Plot(idx2Frame,1)-mean(data2Plot(idx2Frame,1));
           allDataM(i,frames(idx2Frame)) = data2Plot(idx2Frame,2)-mean(data2Plot(idx2Frame,2));
           
           %calculate the std of on particle in one step
           prec(j)  = std(data2Plot(idx2Frame,1));
           %for accuracy we calculate the mean so we can later calculate
           %the difference between steps and thus the accuracy by comparing
           %by the motor step
           acc(j)   = mean(data2Plot(idx2Frame,1));
           precM(j) = std(data2Plot(idx2Frame,2));
           accM(j)  = mean(data2Plot(idx2Frame,2));
           int(j)   = mean(currTrace.intensity(idx2Frame));
           SNR(j)   = mean(currTrace.SNR(idx2Frame));

       end
    end
    plot(currTrace.t,data2Plot(:,1))
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
idx2 =lenTraces>max(lenTraces-period);
tracePlot = tracePlot(idx2);

for i = 1:length(tracePlot)
    currTrace = tracePlot{i};
    currMot = mot{i};
    data2Plot = getData2Plot(currTrace,dim,currMot);
    plot(1:length(data2Plot),data2Plot(:,1))%,5,'filled')
    
end
motPlot = mot{i}*1000;

if ~strcmp(dim,'z')
    motPlot = abs(motPlot);
end

%plot with errorBar
idx = isnan(nanmean(allData,1));
allData(:,idx) =[];
allStd = zeros(nFiles,size(allData,2));
for i = 1:nFiles
    idx = fileRef ==i;
    allStd(i,:) = nanstd(allData(idx,:),1);     
    
end
error = nanmean(allStd,1);

motPlot = motPlot(1:length(error));
motPlot = motPlot - mean(motPlot);
%plot(motPlot,'-r','LineWidth',2);
figure
h= Plotting.shadedErrorBar(1:length(error),motPlot,3*error);

%% function

function [data2Plot] = getData2Plot(currTrace,dim,mot)
    
    switch(dim)
        
        case 'x'
            mot = abs(mot*1000);
            mot = mot-mean(mot);
            currMot = mot(currTrace.t);
            data2Plot(:,1) = currTrace.col-mean(currTrace.col) + mean(currMot);
            data2Plot(:,2) = currTrace.colM-mean(currTrace.colM) + mean(currMot);
            
        case 'y'
            mot = abs(mot*1000);
            mot = mot-mean(mot);
            currMot = mot(currTrace.t);
            data2Plot(:,1) = currTrace.row-mean(currTrace.row) + mean(currMot);
            data2Plot(:,2) = currTrace.rowM-mean(currTrace.rowM) + mean(currMot);
        
        case 'z'
            mot = mot*1000;
            mot = mot - mean(mot);
            currMot = mot(currTrace.t);
            data2Plot(:,1) = currTrace.z-mean(currTrace.z) + mean(currMot);
            data2Plot(:,2) = currTrace.zM-mean(currTrace.zM) + mean(currMot);
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