%%
clear all
clc
close all;
%% get path to zCalibration


path2Cal = 'F:\Data\Leuven Data\2019\06 - June\Roger\29\2D';
calpath ='F:\Data\Leuven Data\2019\06 - June\Roger\29\Laser focus\plane configuration_1';

path = [calpath];

raw.path = path;
raw.ext = '.ome.tif';

%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
info.zMethod    = 'PSFE';%or 'PSFE'

%%
calib = Core.MPMovie(raw,path2Cal,info);
calib.calibrate;
maxFrame = calib.raw.maxFrame(1);

data = calib.getFrame(floor(maxFrame/2));
[Mask] = Misc.getROIFromUser(data(:,:,4));

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

