function makeTraceMovie(data2Store,fullStack,filename,frameRate)

colorMarker = {'w+','y+','g+','c+','b+','r+'};
colorLine = {'w-','y-','g-','c-','b-','r-'};
nFrames = size(fullStack,3);
Fig = figure;

for i = 1:nFrames
    
    cFrame = fullStack(:,:,i);

    imagesc(cFrame)
    hold on
    colormap('gray')
    axis image;
    
    for j = 1 : size(data2Store,3)
        plot(data2Store(i,1,j),data2Store(i,2,j),colorMarker{j})
        plot(data2Store(1:i,1,j),data2Store(1:i,2,j),colorLine{j})
    end
    
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