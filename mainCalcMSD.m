clc ;
clear ;
close all;

%% USER INPUT
expTime = 0.01; %in sec
T = 296; %temperature in Kelvin
R = 0.286; %Radius of particle in um;
fitRDiff = 4; %in nuber of data
minSize = 100; %frames
ext = '.mat';
path = 'D:\Documents\Unif\PhD\2019-Data\Viscosity\500nm-2';
%% Loading
folder = dir(path);
idx = contains({folder.name},'trackResults.mat');
folder(~idx) = [];

f2Load = [folder(1).folder filesep folder(1).name];

tmpData = load(f2Load);
name = fieldnames(tmpData);
data = tmpData.(name{1});

%% Processing
currMov =  data(1).traces;
allHeight = cellfun(@height,currMov(:,1));
idx = allHeight>minSize;
currMov = currMov(idx,1);
allRes = struct('msdX',0,'msdY',0,'msdZ',0,'msdR',0,'tau',0,'DX',0,'DY',0,'DZ',0,'DR',0,...
    'nX',0,'nY',0,'nZ',0,'nR',0);
allRes(length(currMov)).msdX = [];
maxLength = max(allHeight);
allMSDX = zeros(length(currMov),maxLength-1);
allMSDY = allMSDX;
allMSDZ = allMSDY;

for j = 1:length(currMov)
    currPart = currMov{j};

    coord = [currPart.col, currPart.row, currPart.z];
    CM = mean(coord,1);
    coord = coord-CM;

    %in X
    msdx = MSD.calc(coord(:,1)/10^3,'3D');%convert to um;
    allMSDX(j,1:length(msdx)) = msdx;
    DX   = MSD.getDiffCoeff(msdx,fitRDiff,'1D');
    DX   = DX/expTime;%convert from frame to s-1
    nX    = MSD.getViscosity(DX,R,T);

    %inY
    msdy = MSD.calc(coord(:,2)/10^3,'3D');%convert to um;
    allMSDY(j,1:length(msdy)) = msdy;
    DY   = MSD.getDiffCoeff(msdy,fitRDiff,'1D');
    DY   = DY/expTime;%convert from frame to s-1
    nY   = MSD.getViscosity(DY,R,T);


    %inZ
    msdz = MSD.calc(coord(:,3)/10^3,'3D');%convert to um;
    allMSDZ(j,1:length(msdz)) = msdz;
    DZ   = MSD.getDiffCoeff(msdz,fitRDiff,'1D');
    DZ   = DZ/expTime;%convert from frame to s-1
    nZ   = MSD.getViscosity(DZ,R,T);


    %inR
    msdr = MSD.calc(coord/10^3,'3D');%convert to um;
    DR   = MSD.getDiffCoeff(msdr,fitRDiff,'3D');
    DR   = DR/expTime;%convert from frame to s-1
    nR   = MSD.getViscosity(DR,R,T);

    allRes(j).msdX = msdx;% in um^2
    allRes(j).msdY = msdy;
    allRes(j).msdZ = msdz;
    allRes(j).msdR = msdr;
    allRes(j).tau = (1:length(msdx))'*expTime; % in sec


    allRes(j).DX   = DX;% in um^2 /sec
    allRes(j).DY   = DY;% in um^2 /sec
    allRes(j).DZ   = DZ;% in um^2 /sec
    allRes(j).DR   = DR;% in um^2 /sec

    allRes(j).nX   = nX;
    allRes(j).nY   = nY;
    allRes(j).nZ   = nZ;
    allRes(j).nR   = nR;
    allRes(j).num  = length(msdx);
end


filename = [path filesep 'msdRes.mat'];
save(filename,'allRes');
h = msgbox('Data succesfully saved');

