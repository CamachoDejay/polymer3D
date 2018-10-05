expTime = 0.2;
%% 3D tracking Comparison

trace1 = allTraces5K025{1,1};

trace2 = allTraces1K025{188,1};


xAx = [-250 250];
yAx = [xAx];
zAx = [xAx];


figure()

colPlot = trace1.col - mean(trace1.col);
rowPlot = trace1.row - mean(trace1.row);
zPlot   = trace1.z   - mean(trace1.z);
tPlot   = (1:length(trace1.col))*expTime;
patch([colPlot nan(size(colPlot))],[rowPlot nan(size(colPlot))],...
                [zPlot nan(size(colPlot))],[tPlot(:) nan(size(colPlot))],...
                'EdgeColor','interp','FaceColor','none')
view(3);
xlim(xAx);
ylim(yAx);
zlim(zAx);

xlabel('x position (nm)')
ylabel('y position (nm)')
zlabel('z position (nm)')

figure()

colPlot = trace2.col - mean(trace2.col);
rowPlot = trace2.row - mean(trace2.row);
zPlot   = trace2.z   - mean(trace2.z);
tPlot   = (1:length(trace2.col))*expTime;
patch([colPlot nan(size(colPlot))],[rowPlot nan(size(colPlot))],...
                [zPlot nan(size(colPlot))],[tPlot(:) nan(size(colPlot))],...
                'EdgeColor','interp','FaceColor','none')
view(3);
xlim(xAx);
ylim(yAx);
zlim(zAx);


xlabel('x position (nm)')
ylabel('y position (nm)')
zlabel('z position (nm)')
%%

allData{1} = MSDmat1K025;
%allData{2} = MSDmat1K05;
allData{2} = MSDmat5K025;
%allData{4} = MSDmat5K05;
totFrame = length(MSDmat5K05);
timeFrame = (1:totFrame) *expTime;
figure
hold on
color = {'-r','-b','--r','--b'};
for i = 1 : length(allData)
    allMSD = allData{i};
    
    test = allMSD~=0;
    test = sum(test,1);
    [~,idx] = max(test);
    idx2Use = test>=test(idx)*0.9;
    allMSD = allMSD(:,idx2Use);
    allMSD = allMSD(allMSD(:,1)~=0,:);
    allMSD(allMSD==0) = nan;
    
    err = nanstd(allMSD,[],2);
    
    x = timeFrame(1:length(allMSD));
    y = mean(allMSD,2);
    %plot(x,y)
    Plotting.shadedErrorBar(x(:),y(:),err(:),'lineprops',color{i},'logScale',true);
    
end

xlim([0 100])
ylim([0 1100000])
xlabel('Time (s)');
ylabel('RMS displacement nm3')
box('on')