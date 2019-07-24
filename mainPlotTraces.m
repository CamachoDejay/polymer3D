%% 
% This code assumes that some trackRes has been loaded onto matlab
% no clearing
clc
close all

%% User input

minSize = 200;
expTime = 0.01; %sec

%%

figure
hold on
for i = 1:size(trackRes.traces,1)
    
    currTrace = trackRes.traces{i,1};
    
    if height(currTrace) >= minSize
        colPlot = currTrace.col;
        rowPlot = currTrace.row;
        zPlot   = currTrace.z;
        tPlot   = currTrace.t*expTime;
        %plot with time color coding
        patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')

    end        
end

%% plot only a single traced
figure
idx = 105;

currTrace = trackRes.traces{idx,1};
    
colPlot = currTrace.col;
rowPlot = currTrace.row;
zPlot   = currTrace.z;
tPlot   = currTrace.t*expTime;

patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')

view(3)
colorbar
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image