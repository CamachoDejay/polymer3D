%%%%%%%%%%%%%%%%%%%%%%%%% MAIN TRACKING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
% The aim this code is to plot/view tracking data in a nice way. Coding it
% outside of the object allows for a bit of flexibility.

%% USER INPUT

filePath = 'E:\Data\Leuven Data\2019\04 - April\3\XYZ - CS\X';
dim = 'x';
period = 19;
idx2Plot = 1;
%% LOADING

fileName = [filePath filesep 'trackResults.mat'];

trackData = load(fileName);
name = fieldnames(trackData);
trackData = trackData.(name{1});

%% Process Data

traces = trackData.traces(:,1);
[step,mot] = getMotor(trackData,dim);
fileRef = trackData.traces(:,2);
nTraces = length(traces);
lenTraces = cellfun(@height,traces);

maxLen = max(lenTraces);

idx2Mot = [0 period-1:period:maxLen];

precPerTrace = zeros(nTraces,1);
accPerTrace = zeros(nTraces,1);
counter = 0;
precPerStep =[];
accPerStep =[];
motPlot = zeros(maxLen,1);
figure
hold on
for i =1:nTraces
    
    currTrace = traces{i};
    currMot   = mot{i}(2:end);
    currStep  = step{i};
    stepSize = max(currStep)*1000;
    lenCTrace = height(currTrace);
    
    if lenCTrace == maxLen
        
       CM = mean([currTrace.row,currTrace.col,currTrace.z],1);
       
       data2Plot = getData2Plot(currTrace,dim);
       prec = zeros(length(idx2Mot)-1,1);
       acc = prec;
       
       for j = 1: length(idx2Mot)-1
           
           idx = idx2Mot(j)+1:idx2Mot(j+1);
           prec(j) = std(data2Plot(idx));
           acc(j)  = mean(data2Plot(idx));
           
           
       end
       plot(data2Plot)
       accPerStep = [accPerStep; abs(diff(acc))-stepSize];
       precPerStep  = [precPerStep; prec];
       accPerTrace(i) = mean(abs(diff(acc)));
       precPerTrace(i)  = mean(prec);
       counter = counter+1;
    
    end
end
accPerTrace(accPerTrace==0) = [];
precPerTrace(precPerTrace==0)   = [];

precision = mean(abs(accPerTrace-stepSize));
accuracy  = mean(precPerTrace);

fprintf('The average tracking precision is %d nm based on %d traces\n',round(accuracy),counter);
fprintf('The average tracking accuracy is %d nm  based on %d traces \n',round(precision),counter);



figure
subplot(1,2,1)
histogram(precPerStep)
title('Precision');
subplot(1,2,2)
histogram(accPerStep)
title('Accuracy')




%% Plotting
figure 
hold on
%select the traces to plot
idx = cell2mat(fileRef)==idx2Plot;
tracePlot = traces(idx);
lenTraces = lenTraces(idx);
idx2 =lenTraces==max(lenTraces);
tracePlot = tracePlot(idx2);

for i = 1:length(tracePlot)
    currTrace = tracePlot{i};
    currMot = mot{i};
    data2Plot = getData2Plot(currTrace,dim);
    scatter(1:length(data2Plot),data2Plot,'filled')
    
end
motPlot = abs(mot{i})*1000;
motPlot = motPlot(currTrace.frame);
motPlot = motPlot - mean(motPlot);
plot(motPlot,'-r','LineWidth',2);
%% function

function [data2Plot] = getData2Plot(currTrace,dim)

    switch(dim)
        
        case 'x'
            data2Plot = currTrace.col-mean(currTrace.col);
            
        case 'y'
            data2Plot = currTrace.row-mean(currTrace.row);
        
        case 'z'
            
            data2Plot = currTrace.z-mean(currTrace.z);
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