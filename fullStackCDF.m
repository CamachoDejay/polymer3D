% this script loads data of pore sizes and then calculates CDF for a full
% stack
clear
% close all
clc

mainPath = uigetdir;
assert(ischar(mainPath),'User canceled the selection of file, excecution aborted');
%Extract the part of the folder that is a tif file
Folder_Content = dir(mainPath);
idx   = contains({Folder_Content.name},'FullResults.mat');
files2load = Folder_Content(idx);

%%
idx = strfind(mainPath,filesep);
Fname = mainPath(idx(end)+1:end);

nIter = length(files2load);

totalV = [];
for i = 1:nIter

path2data = strcat(files2load(i).folder,filesep);
p2file      = strcat(path2data,files2load(i).name);

load(p2file)

areaV = cat(2,Results.Area);
areaV = areaV(:);

totalV = cat(1,totalV,areaV);
end
[CDF,CCDF] = Misc.getCDF(totalV);

figure()
plot(CCDF.x,CCDF.y)
a = gca;
a.XScale = 'log';
a.YScale = 'log';
title('CCDF for all files in folder')
xlabel('Pore size')
ylabel('CCDF - prob [0 1]')
a.FontSize = 14;
title(Fname)

save([mainPath filesep Fname 'CCDF.mat'],'CCDF') 
return
%%
clear 
close all
clc
load('J:\Project Z-stacks\2K- 1 mg per ml\2K- 1 mg per mlCCDF.mat')
k2v1 = CCDF
load('J:\Project Z-stacks\2K- 0.75 mg per ml\2K- 0.75 mg per mlCCDF.mat')
k2v075 = CCDF
load('J:\Project Z-stacks\2K- 0.5 mg per ml\2K- 0.5 mg per mlCCDF.mat')
k2v05 = CCDF

figure
plot(k2v05.x,k2v05.y)
hold on
plot(k2v075.x,k2v075.y)
plot(k2v1.x,k2v1.y)
hold off
a=gca;
a.XScale = 'log';
a.YScale = 'log';
a.FontSize = 14;
ylim([1e-6 1])
xlim([1e-1 1e4])
legend('0.5','0.75','1.00')
clear
% close all
clc

mainPath = uigetdir;
assert(ischar(mainPath),'User canceled the selection of file, excecution aborted');
%Extract the part of the folder that is a tif file
Folder_Content = dir(mainPath);
idx   = contains({Folder_Content.name},'FullResults.mat');
files2load = Folder_Content(idx);

%%
idx = strfind(mainPath,filesep);
Fname = mainPath(idx(end)+1:end);

nIter = length(files2load);

totalV = [];
for i = 1:nIter

path2data = strcat(files2load(i).folder,filesep);
p2file      = strcat(path2data,files2load(i).name);

load(p2file)

areaV = cat(2,Results.Area);
areaV = areaV(:);

totalV = cat(1,totalV,areaV);
end
[CDF,CCDF] = Misc.getCDF(totalV);

figure()
plot(CCDF.x,CCDF.y)
a = gca;
a.XScale = 'log';
a.YScale = 'log';
title('CCDF for all files in folder')
xlabel('Pore size')
ylabel('CCDF - prob [0 1]')
a.FontSize = 14;
title(Fname)

save([mainPath filesep Fname 'CCDF.mat'],'CCDF') 
return
%%
clear 
close all
clc
load('J:\Project Z-stacks\2K- 1 mg per ml\2K- 1 mg per mlCCDF.mat')
k2v1 = CCDF
load('J:\Project Z-stacks\2K- 0.75 mg per ml\2K- 0.75 mg per mlCCDF.mat')
k2v075 = CCDF
load('J:\Project Z-stacks\2K- 0.5 mg per ml\2K- 0.5 mg per mlCCDF.mat')
k2v05 = CCDF

figure
plot(k2v05.x,k2v05.y)
hold on
plot(k2v075.x,k2v075.y)
plot(k2v1.x,k2v1.y)
hold off
a=gca;
a.XScale = 'log';
a.YScale = 'log';
a.FontSize = 14;
ylim([1e-6 1])
xlim([1e-1 1e4])
legend('0.5','0.75','1.00')
