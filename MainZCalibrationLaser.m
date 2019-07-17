%%
clear all
clc
close all;
%% get path to zCalibration


path2Cal = 'C:\Users\abaku\Desktop\TRACKING\FEATUREDETECTION\plane configuration_7\2DCal\';
calpath ='C:\Users\abaku\Desktop\TRACKING\FEATUREDETECTION\plane configuration_7\';

path = [path2Cal ];

raw.path = path;
raw.ext = '.ome.tif';

%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
info.zMethod    = 'Intensity';%or 'PSFE'

%%
calib = Core.MPMovie(raw,calpath,info);
calib.calibrate;
maxFrame = calib.raw.maxFrame(1);

data = calib.getFrame(floor(maxFrame/2));
[Mask] = Misc.GetROIFromUserInput(data(:,:,4),calib.cal2D.file.ROI );

chData1las = zeros(size(data,1), size(data,2), size(data,3)/2, maxFrame);
chData2las = zeros(size(data,1), size(data,2), size(data,3)/2, maxFrame);
chData1nor = zeros(size(data,1), size(data,2), size(data,3)/2, maxFrame);
chData2nor = zeros(size(data,1), size(data,2), size(data,3)/2, maxFrame);
[zStep, zPosMotor] =  calib.getZPosMotor ;

for i=1:maxFrame
    data = calib.getFrame(i);

    for j=1:size(data,3)/2
       chData1las(:,:,j,i) =  Mask.*data(:,:,j);
       chData1nor(:,:,j,i) = ~Mask.*data(:,:,j);
    end
    for j=size(data,3)/2+1:size(data,3)
       chData2las(:,:,j-4,i) =  Mask.*data(:,:,j);
       chData2nor(:,:,j-4,i) = ~Mask.*data(:,:,j);
    end
end
[ focus_metLas, in_focusLas, fitLas ] = mpSetup.cali.getFocusMetric( chData1las,chData2las, zPosMotor, zPosMotor );
[ focus_metNor, in_focusNor, fitNor ] = mpSetup.cali.getFocusMetric( chData1nor,chData2nor, zPosMotor, zPosMotor );

distanceToPlanes = [in_focusLas.zpos]- [in_focusNor.zpos];

close all

