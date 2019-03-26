function makeTraceMovie(data2Store,fullStack,filename,frameRate,scaleBar,tail)


color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780  0.1840],[1 1 1]};

nFrames = size(fullStack,3);
if isempty(tail)
    tail = nFrames;
end
Fig = figure;

for i = 1:nFrames %loop over frames
    %plot original image
    cFrame = fullStack(:,:,i);

    imagesc(cFrame)
    hold on
    colormap('gray')
    axis image;
    
    %loop over particles and plot a trace with "tail" number of data point that precede the current point
    for j = 1 : size(data2Store,3) 
        startIdx = i - tail;
        startIdx(startIdx<1) = 1;
        plot(data2Store(startIdx:i,1,j),data2Store(startIdx:i,2,j),'Color',color{j})
        plot(data2Store(i,1,j),data2Store(i,2,j),'+', 'MarkerEdgeColor',color{j});
       
    end
    
    %add scale  bar
    scaleBarPx = scaleBar/95*1000;
    x = size(cFrame,2)-scaleBarPx-(0.05*size(cFrame,2)):size(cFrame,2)-0.05*size(cFrame,2);
    y = ones(1,length(x))*size(cFrame,1)-0.05*size(cFrame,2);
    text(mean(x),mean(y)-0.05*size(cFrame,1),[num2str(scaleBar) ' µm'],'HorizontalAlignment','center','Color','white','fontWeight','bold','fontSize',14);
    plot(x,y,'-w','LineWidth',3);
    
    set(gca,'visible','off');
    set(gcf,'color','w');
    drawnow;

    
    frame = getframe(Fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if i == 1

        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'loopcount',inf);

    else

        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'writemode','append');

    end
    
    hold off
    clf;
end