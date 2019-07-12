% this script is a work in progress, the plan is to fix here the quick
% callibration of the channels. Basically I want to be able to lead the
% movies and built a dataset of [xDim,yDim,channel,frame] where all
% channels represent the same region of the image to a pixel resolution.
clear
close all
clc

file.path  = 'F:\Data\Leuven Data\2019\07 - July\20190708\2D';
file.ext   = '.ome.tif';
info.runMethod = 'load';
info.nChan = 4;
%% 
calib = Core.MPPlaneCalibration(file,info);

calib.retrieveMovies;
calib.calcIndivCal; 
calib.calcCombinedCal;

calib.showCal(1)
calib.offTarget;
calib.save;
%% PREV VERSION
fPath = 'E:\Data\Leuven Data\2018\08-Aug\29\PlaneCalibration_25';
fName = 'PlaneCalibration_25_MMStack_Pos0.ome.tif';

fPath = [fPath filesep fName];

% Calculate calibration
[cal] = mpSetup.cali.calculate(fPath, false);

%TODO: Improve the calculation of distance between the plane, Move it to
%calculate ?
zFocus = cell2mat({cal.inFocus.zpos});
Fit = cal.fit;
% zFocus2 = zeros(1,size(cal.focusMet,2));
% Fit = zeros(size(cal.focusMet));
% for k=1:size(cal.focusMet,2)
%     [out,Fit(:,k)] = SimpleFitting.gauss1D(cal.focusMet(:,k),cal.Zpos);
%     zFocus(k) = out(2);
% end
zFocus = zFocus(cal.neworder); %give the right order to the channels
distBetweenCamPlanes = abs(mean(diff(zFocus(1:2:end))) + mean(diff(zFocus(2:2:end))))/2;
target    = distBetweenCamPlanes/2;
distBetweenPlane = diff(zFocus);
offTarget1 = distBetweenPlane - target;
offTarget = mean(abs(offTarget1));

message = sprintf('The difference between the target and the current plane conformation \nis %d nm',round(offTarget*1000));
disp(message);
%%
 figure()
            hold on
            color = rand(8,3);
            height = max(max(cal.focusMet));
            for i = 1 : size(cal.focusMet,2)
                [~,idx] = max(Fit(:,i));
                scatter(cal.Zpos(:),cal.focusMet(:,i),[],color(i,:),'filled')
                plot(cal.Zpos(:),Fit(:,i),'Color', color(i,:),'LineWidth',2.5)
                
                
                y = 1:height;
                x = ones(1,length(y))*zFocus(i);
                plot(x(:),y(:),'k--');

            end
            ylim([min(min(cal.focusMet)), max(max(cal.focusMet))]);
            xlim([min(cal.Zpos), max(cal.Zpos)]);
            title('Setup Plane Calibration');
            
            hold off

%%
% load and calibrate, when applied to the calibration data then we should
% be able to demonstrate that it works
[data] = mpSetup.loadAndCal( fPath, cal );

% figure to show that it works
f = round(mean(cat(1,cal.inFocus.frame)));
xl = [97,206];
yl = [116,229];
for i =1:8
    subplot(2,4,i)
    imagesc(data(:,:,i,f))
    axis image
%     xlim(xl)
%     ylim(yl)
%     hold on
    grid on
    
    title(['Ordered channel ' num2str(i)])
    a = gca;
    a.XTickLabel = [];
    a.YTickLabel = [];
end



