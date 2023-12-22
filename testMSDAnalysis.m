% This code aim at analyzing trajectories and extracting useful metric out
% of the traces (e.g. Diffusion coefficient, viscosity,...)
%
% This code works well with the trackRes output that most of the tracking
% codes output, it just need to be given the path to the folder where the
% trackRes.mat is store, and a few input which are summarized in the frame
% below, the rest is handled by the code internally.
%
clc ;
clear ;
close all;
%% USER INPUT
expTime = 0.02; %in sec
T = 293; %temperature in Kelvin
R = 0.075; %Radius of particle in um;
fitRDiff = 4; %in Fraction of the data
fitRConf = 0.1;%in Fraction of the data
minSize = 10;
ext = '.mat';
path = 'E:\Results\SPT trapping\Trapping 640\pH 3';

%% Loading
folder = dir(path);
idx = contains({folder.name},'trackRes.mat');
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

for i = 1: size(data,2)
   
    warning('only analyzing 2D now')
    %get CM of all particles
    coord = zeros(length(currMov),2);
    for j = 1:length(currMov)
        
        coord(j,:) = [mean(currMov{j}.col),mean(currMov{j}.row)];%,mean(currMov{j}.z)];
                
    end
    CM = mean(coord,1);
    
    allRes = struct('msd',zeros(length(coord),1),'tau',zeros(length(coord),1),...
        'D',0,'n',0,'alpha',0,'rConf',0,'stiff',0);
    allRes(length(currMov)).msd = zeros(length(coord),1);
    
    for j = 1:length(currMov)
        currPart = currMov{j};
        
        coord = [currPart.col, currPart.row];
        
        coord = coord-CM;
        
        msd = MSD.calc(coord/10^3);%convert to um;
        D   = MSD.getDiffCoeff(msd,1:length(msd),fitRDiff,'2D');
        D   = D/expTime;%convert from frame to s-1
        n   = MSD.getViscosity(D,R,T);
        
        alpha = MSD.getDiffTypeAlpha(msd,expTime);
        if alpha < 1
            rConf = MSD.getConfRad(msd,fitRConf,expTime);
        else
            rConf = NaN;
        end
        stiff = MSD.getTrapStiffness(coord*10^(-9),T);%convert to meter
        
        allRes(j).msd = msd; % in um^2
        allRes(j).tau = (1:length(msd))'*expTime; % in sec
        allRes(j).D   = D;% in um^2 /sec
        allRes(j).n   = n;% in um^2 /sec 
        allRes(j).alpha = alpha;%dimension less
        allRes(j).rConf = rConf;% in um
        allRes(j).stiff = stiff;% in pN/um
    end
    
    allData(i).fileName = data(i).fileName;
    allData(i).msdRes   = allRes;
end

filename = [path filesep 'msdRes.mat'];
save(filename,'allData');
h = msgbox('Data succesfully saved');
