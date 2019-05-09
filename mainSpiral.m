%%%%%%%%%%%%%%%%%%%%%%%%% MAIN TRACKING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
% The aim this code is to plot/view tracking data in a nice way. Coding it%
% outside of the object allows for a bit of flexibility.                  %
clear;
close all;
clc;
%% USER INPUT

filePath = 'E:\Data\Leuven Data\2019\04 - April\3\Spirals';
ignoreF = 6;
dim = 'z';
period = 20;
idx2Plot = 3;
stepApplied = 200;
longTraceThresh = 10;
%% LOADING

fileName = [filePath filesep 'trackResults.mat'];

trackData = load(fileName);
name = fieldnames(trackData);
trackData = trackData.(name{1});

%% Aligning all curves
lengthTraces = cellfun(@height,trackData.traces(:,1));
idx2LTraces = lengthTraces>longTraceThresh;

traces = trackData.traces(idx2LTraces,1);
xMot   = trackData.traces(idx2LTraces,4);
yMot   = trackData.traces(idx2LTraces,6);
zMot   = trackData.traces(idx2LTraces,8);
file   = trackData.traces(idx2LTraces,2);
frame = 1:length(xMot{1});
alignedTraces = cell(size(traces));
Spiral = alignedTraces;
%alignTrace
for i = 1:length(traces)
    cTrace = traces{i};
    cFile  = file{i};
    cTrace.row = -cTrace.row;
    cTrace.col = -cTrace.col;
    
    cXMot = xMot{i}*1000;
    cYMot = yMot{i}*1000;
    cZMot = zMot{i}*1000;
    
    Spiral{i} = [cXMot - mean(cXMot),cYMot - mean(cYMot),cZMot - mean(cZMot),frame(:)];
    
    cFrames = cTrace.frame;
    alignedTraces{i} = [cTrace.col - mean(cTrace.col) + mean(Spiral{i}(cFrames,1)),...
        cTrace.row - mean(cTrace.row) + mean(Spiral{i}(cFrames,2)),...
        cTrace.z - mean(cTrace.z) + mean(Spiral{i}(cFrames,3)),cFrames,ones(length(cFrames),1)*cFile];
      
end
cSpiral = Spiral{1};
%% Processing
%precision preallocation
xPrec = zeros(length(Spiral{1}),1);
yPrec = zeros(length(Spiral{1}),1);
zPrec = zeros(length(Spiral{1}),1);
%Accuracy preallocation
xAcc = zeros(length(Spiral{1}),1);
yAcc = zeros(length(Spiral{1}),1);
zAcc = zeros(length(Spiral{1}),1);

trace2Proc = cell2mat(alignedTraces);
for i = 1:length(cSpiral)
    cTrace2Proc = trace2Proc(trace2Proc(:,4)==i,:);
    
    if ~isempty(cTrace2Proc)
        xPrec(i) = std(cTrace2Proc(:,1));
        yPrec(i) = std(cTrace2Proc(:,2));
        zPrec(i) = std(cTrace2Proc(:,3));
        
        xAcc(i)  = mean(abs(cTrace2Proc(:,1)-cSpiral(i,1)));
        yAcc(i)  = mean(abs(cTrace2Proc(:,2)-cSpiral(i,2)));
        zAcc(i)  = mean(abs(cTrace2Proc(:,3)-cSpiral(i,3)));
        
        
    end
  
end
xPrec(xPrec==0) = [];
yPrec(yPrec==0) = [];
zPrec(zPrec==0) = [];
xAcc(xPrec==0)  = [];
yAcc(xPrec==0)  = [];
zAcc(xPrec==0)  = [];

%% Plotting per file
figure
hold on
fileIdx = cell2mat(file);
idx = fileIdx==idx2Plot;

traces2Plot = alignedTraces(idx);

for i = 1:length(traces2Plot)
    cTrace = traces2Plot{i};
    scatter3(cTrace(:,1),cTrace(:,2),cTrace(:,3),10,cTrace(:,3),'filled');
    plot3(cTrace(:,1),cTrace(:,2),cTrace(:,3))
end

%% Plotting
%initialize sphere
nfacets = 15;
res = 50;
[sx,sy,sz]= sphere(nfacets);
figure
hold on
plot3(cSpiral(:,1),cSpiral(:,2),cSpiral(:,3),'LineWidth',3);

for i = 1:length(cSpiral(:,1))
    Sx = sx*res+cSpiral(i,1);
    Sy = sy*res+cSpiral(i,2);
    Sz = sz*res+cSpiral(i,3);
    surf(Sx,Sy,Sz)
    
end
hold on
for i = 1: length(traces)

    cTrace = alignedTraces{i};  
    scatter3(cTrace(:,1),cTrace(:,2),cTrace(:,3),10,cTrace(:,3),'filled');
    
end
colormap('jet')
hold off

fprintf('The mean precision in\n x =%d, y=%d, z=%d\n',median(xPrec),median(yPrec),median(zPrec));
fprintf('The mean accuracy in\n x =%d, y=%d, z=%d\n',median(xAcc),median(yAcc),median(zAcc));
