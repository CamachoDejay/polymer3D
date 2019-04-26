%%%%%%%%%%%%%%%%%%%%%%%%% MAIN TRACKING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
% The aim this code is to plot/view tracking data in a nice way. Coding it
% outside of the object allows for a bit of flexibility.

%% USER INPUT

filePath = 'E:\Data\Leuven Data\2019\04 - April\3\XYZ - CS\X';
dim = 'x';
period = 19;

%% LOADING

fileName = [filePath filesep 'trackResults.mat'];

trackData = load(fileName);
name = fieldnames(trackData);
trackData = trackData.(name{1});

%% Process Data

traces = trackData.traces(:,1);
[step,mot] = getMotor(trackData,dim);

nTraces = length(traces);
lenTraces = cellfun(@height,traces);

maxLen = max(lenTraces);

idx2Mot = [0 period-1:period:maxLen];

accPerTrace = zeros(nTraces,1);
precPerTrace = zeros(nTraces,1);
counter = 0;
accPerStep =[];
precPerStep =[];
figure
hold on
for i =1:nTraces
    
    currTrace = traces{i};
    currMot   = mot{i};
    currStep  = step{i};
    stepSize = max(currStep)*1000;
    lenCTrace = height(currTrace);
    
    if lenCTrace == maxLen
        
       CM = mean([currTrace.row,currTrace.col,currTrace.z],1);
       
       data2Plot = getData2Plot(currTrace,dim);
       acc = zeros(length(idx2Mot)-1,1);
       prec = acc;
       
       for j = 1: length(idx2Mot)-1
           
           idx = idx2Mot(j)+1:idx2Mot(j+1);
           acc(j) = std(data2Plot(idx));
           prec(j)  = mean(data2Plot(idx));          
           
       end
       plot(data2Plot)
       precPerStep = [precPerStep; abs(diff(prec))-stepSize];
       accPerStep  = [accPerStep; acc];
       precPerTrace(i) = mean(abs(diff(prec)));
       accPerTrace(i)  = mean(acc);
       counter = counter+1;
    
    end
end
precPerTrace(precPerTrace==0) = [];
accPerTrace(accPerTrace==0)   = [];

precision = mean(abs(precPerTrace-stepSize));
accuracy  = mean(accPerTrace);

fprintf('The average tracking accuracy is %d nm based on %d traces\n',round(accuracy),counter);
fprintf('The average precision is %d nm  based on %d traces \n',round(precision),counter);

%% Plotting



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