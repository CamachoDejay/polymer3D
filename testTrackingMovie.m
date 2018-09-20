%%
clc 
close all
clear 

path2ZCal = 'E:\Data\Leuven Data\2018\06-June\27\ZCal - NormObjCorr';
path2SRCal = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - normObjCorrPSFE';

path2File = 'E:\Data\Leuven Data\2018\06-June\27\ZCal - NormObjCorr\zStackFluoBeads200_PIC1_270618__5';
path2Cal = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - normObjCorrPSFE\zStackFluoBeads200_S3_270618__1';

detectParam.delta = 6;
detectParam.chi2 = 80;
%%
calib = Core.MPCalibration(path2Cal);

%%

MPTrackMov = Core.MPTrackingMovie(path2File,calib.getCal,path2SRCal,path2ZCal);

%% Detection

MPTrackMov.giveInfo
%find candidate
MPTrackMov.findCandidatePos(detectParam);

%fit position
MPTrackMov.SRLocalizeCandidate;

%% Data correction
rot = true;
refPlane = 5;
MPTrackMov.applySRCal(rot,refPlane);
%% e-Z transformation
MPTrackMov.applyZCal;

%% Plane consolidation

MPTrackMov.consolidatePlanes
%% Super resolve
val2Use = 'bestFocus';
MPTrackMov.superResolve(val2Use);
%% plot
frames = 1:100;

MPTrackMov.showCorrLoc();

%% showFrame

MPTrackMov.showFrame(80);
%MPTrackMov.showParticle;

%% tracking
trackParam.euDistXY = 250;
trackParam.euDistZ  = 300;
MPTrackMov.trackParticle(trackParam);
%% plot
MPTrackMov.showTraces;

%% eval accuracy
MPTrackMov.evalAccuracy

%% Susana's figure Movies
raw = MPTrackMov.getRaw;
traces = MPTrackMov.getTraces;
%roiRadius = 10;
currentTraces = traces {1};
frameRate = 5;
mainPos = [round(mean(currentTraces.row)/95) round(mean(currentTraces.col(1)/95))];
frames = 250;
for i = 1:8
    currentPlane = MPTrackMov.getPlane(i);
    ROI = currentPlane;
   % ROI = currentPlane(mainPos(1)-roiRadius:mainPos(1)+roiRadius,...
      %  mainPos(2)-roiRadius:mainPos(2)+roiRadius,:);
    mov = struct('cdata', cell(1,frames), 'colormap', cell(1,frames));
    Fig = figure;
  %to get as less white border as possible
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    for j = 1:frames
        
        gcf;
        
        imagesc(ROI(:,:,j))
        hold on
        %scale bar
        x = size(ROI,2)-90:size(ROI,2)-50;
        y = ones(1,length(x))*size(ROI,1)-40;
        plot(x,y,'-w','LineWidth',5);
        caxis([min(min(min(ROI))) max(max(max(ROI)))]);
        axis image;
        set(gca,'visible','off');
        colormap('hot')
        drawnow;
        
        hold off
        
        mov(j) = getframe(Fig);

    end
    ext='.mp4';
    filename=sprintf('%s%splane%d%s', path2File,'\',i,ext);
    v = VideoWriter(filename,'MPEG-4');
    v.FrameRate = frameRate;
    open(v)
    writeVideo(v,mov);
    close(v)
end
%% 
sizeMarker = 5;
 Fig = figure;
 
 xAx = [-500, 500];
 yAx = [-500, 500];
 zAx = [-1000, 1000];
 mov = struct('cdata', cell(1,size(currentTraces.row,1)), 'colormap', cell(1,size(currentTraces.row,1)));
 hold on
 
    for j = 1:size(currentTraces.row,1)
        col = currentTraces.col(j) - mean(currentTraces.col);
        row = currentTraces.row(j) - mean(currentTraces.row);
        z = currentTraces.z(j) - mean(currentTraces.z);
        scatter3(col,row,z,sizeMarker,z,'filled')
        xlim(xAx)
        ylim(yAx)
        zlim(zAx)
        
        view(3);
        %colormap('hot')
        mov(j) = getframe(Fig);
        
        xlabel('x Position (nm)');
        ylabel('y Position(nm)');
        zlabel('z Position(nm)');
    end

 ext='.mp4';
    filename=sprintf('%s%sTrackingNorm%d%s', path2File,'\',i,ext);
    v = VideoWriter(filename,'MPEG-4');
    v.FrameRate = frameRate;
    open(v)
    writeVideo(v,mov);
    close(v)