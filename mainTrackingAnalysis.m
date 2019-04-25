%%%%%%%%%%%%%%%%%%%%%%%%% MAIN TRACKING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
% The aim this code is to plot/view tracking data in a nice way. Coding it
% outside of the object allows for a bit of flexibility.

%% USER INPUT

filePath = 'E:\Data\Leuven Data\2019\04 - April\3\XYZ - CS\X';
dim = 'x';
period = 19;
stepSize = 200;
%% LOADING

fileName = [filePath filesep 'trackResults.mat'];

trackData = load(fileName);
name = fieldnames(trackData);
trackData = trackData.(name{1});

%% Process Data

traces = trackData.traces;
nTraces = length(traces);
lenTraces = cellfun(@height,traces);

maxLen = max(lenTraces);

idx2Mot = [0 period-1:period:maxLen];

accPerTrace = zeros(nTraces,1);
precPerTrace = zeros(nTraces,1);

for i =1:nTraces
    
    currTrace = traces{i};
    lenCTrace = height(currTrace);
    
    if lenCTrace == maxLen
        
       CM = mean([currTrace.row,currTrace.col,currTrace.z],1);
       
       data2Plot = getData2Plot(currTrace,CM,dim);
       acc = zeros(length(idx2Mot)-1,1);
       prec = acc;
       
       for j = 1: length(idx2Mot)-1
           
           idx = idx2Mot(j)+1:idx2Mot(j+1);
           acc(j) = std(data2Plot(idx));
           prec(j)  = mean(data2Plot(idx));          
           
       end
       
       precPerTrace(i) = mean(abs(diff(prec)));
       accPerTrace(i)  = mean(acc);
       
    
    end
end
precPerTrace(precPerTrace==0) = [];
accPerTrace(accPerTrace==0)   = [];

precision = mean(abs(precPerTrace-stepSize));
accuracy  = mean(accPerTrace);

fprintf('The average tracking accuracy is %d nm \n',round(accuracy));
fprintf('The average precision is %d nm \n',round(precision));

function [data2Plot] = getData2Plot(currTrace,CM,dim)

    switch(dim)
        
        case 'x'
            data2Plot = currTrace.col-CM(2);
            
        case 'y'
            data2Plot = currTrace.row-CM(1);
        
        case 'z'
            
            data2Plot = currTrace.z-CM(3);
    end

end
    