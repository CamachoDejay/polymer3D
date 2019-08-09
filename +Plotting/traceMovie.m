%% User Input

path2Save = 'D:\Documents\Unif\PhD\Papers\04 - Particle Tracking\Figure';

ext = '.gif';
filename=sprintf('%s%s3Dejection%s', path2Save,filesep,ext);

data = allTraces;
Trailing = 100;
frameRate = 50;
yLimit = [2000 10000];
xLimit = [12500 22500];
zLimit = [-4500 -0];
%% test data type

if istable(data{1})
    dataType = 'table';
elseif ismatrix(data{1})
    dataType = 'matrix';
end
%%
%function traceMovie(data,frameRate)
%Find maxima and clean data
maxFr = zeros(length(data),1);
idx2Delete = zeros(length(data),1);
for i = 1:length(data)
    currData = data{i};
    switch dataType
        case 'table'
            maxFr(i) = max(currData.t);
            testY = and(currData.row > yLimit(1), currData.row < yLimit(2));
            testX = and(currData.col > xLimit(1), currData.col < xLimit(2));
        case 'matrix'
            maxFr(i) = max(currData(:,4));
            testY = and(currData(:,2) > yLimit(1), currData(:,2) < yLimit(2));
            testX = and(currData(:,1) > xLimit(1), currData(:,1) < xLimit(2));
    end
    test = testX.*testY;

    if all(test==0)
        idx2Delete(i) = i;
    end
end
idx2Delete(idx2Delete==0) = [];

data(idx2Delete) = [];

%% Localization Density image

ZStep = 500;
XBin = 500;
YBin = 500;
countBound = [0 100];

%Put all the data into a single table
tableData = data{1};
for i = 2:length(data)
    tableData = [tableData; data{i}];

end

zLim = [min(tableData.z) max(tableData.z)];
zBin = zLim(1):ZStep:zLim(2);

%Binning along Z
zBinnedData = cell(length(zBin)-1,1);
for i=1:length(zBin)-1
    
    zBinnedData{i} = tableData(and(tableData.z>zBin(i),tableData.z<zBin(i+1)),:);
end

xLim = [min(tableData.col) max(tableData.col)];
xBin = xLim(1):XBin:xLim(2);

yLim = [min(tableData.row) max(tableData.row)];
yBin = yLim(1):YBin:yLim(2);

histData = zeros(numel(xBin)-1,numel(yBin)-1,length(zBinnedData));
for i = 1:length(zBinnedData)
    currData = zBinnedData{i};
     h = histogram2(currData.col,currData.row,xBin,yBin);
    histData(:,:,i) = h.BinCounts;
end
figure
hold on
for i = 1:size(histData,3)
   subplot(2,round(size(histData,3)/2),i)
   imagesc((yBin(1:end-1)-yBin(1))/1000,(xBin(1:end-1)-xBin(1))/1000,histData(:,:,i));
   caxis(countBound);
   title( ['Z = ' num2str(round(zBin(i)))])
   axis image; 
   colormap('jet')
end

%% Make Movie
 

    %get max Frame
    maxFr = max(maxFr);
%     [X,Y,Z] = sphere(15);
%        SX = X*20+data2Plot.col(currData.frame==i);
%                 SY = Y*20+data2Plot.row(currData.frame==i);
%                 SZ = Z*20+data2Plot.z(currData.frame==i);
    Fig = figure;
   view(3);
     xlim(xLimit);
     ylim(yLimit);
    xlim manual;
    ylim manual;
    gcf;
    hold on
    
    for i = 1:200
       
        for j = 1:length(data)
            currData = data{j};
            switch dataType
                case 'table'
                    idx2Frame = currData.frame==i;
                    idx = i-100:i;
                    idx = ismember(currData.frame,idx);
                case 'matrix'
                    
                    idx2Frame = currData(:,4) == i;
                    idx = i-100:i;
                    idx = ismember(currData(:,4),idx);
                    
            end
            if ~all(idx==0)
            
                data2Plot = currData(idx,:);
                xlim(xLimit);
                ylim(yLimit);
                zlim(zLimit)
                xlim manual;
                ylim manual;
                zlim manual;
                view(3);
                gcf;
                hold on
                switch dataType
                    case 'table'
                        plot3(data2Plot.col,data2Plot.row,data2Plot.z,'color',[0 0 1])             
                        scatter3(currData.col(idx2Frame),currData.row(idx2Frame),currData.z(idx2Frame),40,[1 0 0],'filled')
                    case 'matrix'
                        plot3(data2Plot(:,1),data2Plot(:,2),data2Plot(:,3),'color',[0 0 1])             
                        scatter3(currData(idx2Frame,1),currData(idx2Frame,2),currData(idx2Frame,3),40,[1 0 0],'filled')
                   
                end
            end
        end
        drawnow;
        frame = getframe(Fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);

        if i == 1

            imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'loopcount',inf);

        else

            imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'writemode','append');

        end
        clf;
    end

    
    %% make Last frame
    figure
    hold on
    for i = 1:length(data)
        
         currData = data{i};
         plot3(currData.col,currData.row,currData.z) 
         
    end

    drawnow;
    frame = getframe(Fig);
    im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'writemode','append');
