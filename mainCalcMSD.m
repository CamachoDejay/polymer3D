clc ;
clear ;
close all;

%% USER INPUT
expTime = 0.01; %in sec
T = 293; %temperature in Kelvin
R = 0.125; %Radius of particle in um;
fitRDiff = 0.05; %in Fraction of the data
fitRConf = 0.1;%in Fraction of the data
ext = '.mat';
path = 'D:\Documents\Unif\PhD\2019-Data\09 - Sep\2P';
%% Loading
folder = dir(path);
idx = contains({folder.name},'trackRes.mat');
folder(~idx) = [];

f2Load = [folder(1).folder filesep folder(1).name];

tmpData = load(f2Load);
name = fieldnames(tmpData);
data = tmpData.(name{1});

%% Processing

for i = 1: size(data,2)
    currMov =  data(i).traces;
      
    allRes = struct('msd',zeros(length(coord),1),'tau',zeros(length(coord),1),...
        'D',zeros(length(coord),1),'n',zeros(length(coord),1));
    allRes(length(currMov)).msd = zeros(length(coord),1);
     
    for j = 1:length(currMov)
        currPart = currMov{j};
        
        coord = [currPart.col, currPart.row, currPart.z];
        CM = mean(coord,1);
        coord = coord-CM;
        
        msd = MSD.calc(coord/10^3,'3D');%convert to um;
        D   = MSD.getDiffCoeff(msd,fitRDiff,'3D');
        D   = D/expTime;%convert from frame to s-1
        n   = MSD.getViscosity(D,R);
        
        allRes(j).msd = msd; % in um^2
        allRes(j).tau = (1:length(msd))'*expTime; % in sec
        allRes(j).D   = D;% in um^2 /sec
        allRes(j).n   = n;% in um^2 /sec
        
     end
    
    allData(i).fileName = data(i).fileName;
    allData(i).msdRes   = allRes;
end

filename = [path filesep 'msdRes.mat'];
save(filename,'allData');
h = msgbox('Data succesfully saved');

