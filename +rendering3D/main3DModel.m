% The aim of this code si to generate a 3D view of the sample making an Iso
% surface of the boundary between 0 and 1 in the binary data. After the
% isosurface is made, the isosurface is smoothen to improve the
% visualization
%
% !!!!!!! the code require to install a matlab compiler for c/c++ that can
% be found at:
% https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler
%
% The smoothing of the data was NOT WRITTEN by the authors, the code written
% in c is already included in the folder but more information can be found
% at: %https://www.mathworks.com/matlabcentral/fileexchange/26710-smooth-triangulated-mesh

% HOW TO USE:
% The code is run, it will open the file selector from window, the user can
% then choose tif files of their binary images (segmented stacks). Only
% binary images are supported. The code will then calculate the isosurface
% and the smoothing.Optionnally, the code will represent sphere where pores
% were detected and make an additional plot with the pores and their 
% connectivity. Finally, the code will display the result and save the figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boris Louis (https://github.com/BorisLouis)                             %
% Website : Boris Louis: https://borislouis.github.io/                    %
% Susana Rocha (https://github.com/SusanaRocha)                           %
% Website: https://susanarocha.github.io/                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
clc;
%compile c code for smoothing
mex +rendering3D\smoothpatch_curvature_double.c -v
mex +rendering3D\smoothpatch_inversedistance_double.c -v
mex +rendering3D\vertex_neighbours_double.c -v
%% User Input
smoothMode = 0; % 0 gives smoother edges, more details in smoothpatch.m
nIterations = 10; %number of smoothing iterations
pxSizeXY = 180; %in nm
pxSizeZ = 404; %in nm
%colors of the plot
colorModel = [0.6,0,0];%replace by : 'Z' for zcoloring
colorPores = [0.4 0.8 0.8];

%Select area to plot
xRange = [30 196];%px
yRange = [30 196];%px

%more parameter are available and can be checked within the smoothpatch
%function, the default parameter put in this implementation were satisfying
%to us.
%% load data
%open file selector
[name,path] = uigetfile('*.tif','Please choose a binary file to model');
dataPoreSize = [path filesep 'PoreSize-Results'];
%get the path to data
path2Data  = [path filesep name];
%extract filename without extension
[~,namenoExt,~] = fileparts(path2Data);

warning('off')
fInfo = Load.Movie.tif.getinfo(path2Data);
%loading the data
data2Render = Load.Movie.tif.getframes(path2Data,1:fInfo.Frame_n);
warning('on');
% check if the data is indeed binary
test = unique(data2Render);
assert(length(test)<=2,'The input data is not binary, please choose a binary file e.g. segmented data');
%if exist get the poreData
if isfolder(dataPoreSize)
    folderContent = dir(dataPoreSize);
    idx = contains({folderContent.name},'.mat');
    file2Load = [folderContent(idx).folder filesep folderContent(idx).name];
    if isfile(file2Load)
        load(file2Load);
        pores=true;
    else
        pores = false;
        
    end
else
    pores = false;
end

%% Preparing data
%from the analysis file
data2Render = data2Render(yRange(1):yRange(2),xRange(1):xRange(2),:);
%convert range to um
xRange = xRange * pxSizeXY*1e-3;
yRange = yRange * pxSizeXY*1e-3;

if pores
    % extract the data within the user input border
    poresCoord=allData(1).pores3D.ctr_ext;
    idxX = and(poresCoord(:,1)>=xRange(1),poresCoord(:,1)<=xRange(2));
    idxY = and(poresCoord(:,2)>=yRange(1),poresCoord(:,2)<=yRange(2));
    idx = logical(idxX.*idxY);
    poresCoord = poresCoord(idx,:);
    %extract the radius accordingly
    extRad=allData(1).pores3D.extRad;
    extRad = extRad(idx);
    
    %extract connectivity
    connID = allData(1).pores3D.connID;   
    connID = connID(idx);
    idxF = find(idx==1);
    %clean connectivity from out-of-range pores that are connected to
    %in-range ones.
    for i = 1 : length(connID)
        list = connID{i};
        idx = ~ismember(list,idxF);
        list(idx) = [];
        connID{i} = list;
        
    end
else
    warning('no pores data were found so they will not be plotted please run mainPore Size calc of the file of interest first.');
end
%Here we decided to take the full Z axis
zRange = [0 max(poresCoord(:,3))];
%% Making isosurface
iSurface = isosurface(data2Render,1/2);
% smoothing using compiled c code
smoothISurface = rendering3D.smoothpatch(iSurface,0,nIterations);
%comnvert to px
smoothISurface.vertices(:,1) = (smoothISurface.vertices(:,1))*pxSizeXY*1e-3;
smoothISurface.vertices(:,2) = (smoothISurface.vertices(:,2))*pxSizeXY*1e-3;
smoothISurface.vertices(:,3) = (smoothISurface.vertices(:,3))*pxSizeZ *1e-3;

%% Displaying network model
%z-coloring
if strcmpi(colorModel,'Z')
    colorModel = smoothISurface.vertices(:,3)/max(smoothISurface.vertices(:,3));
    zColor = true;
else
    zColor = false;
end
%Plot the network with Z coloring or unique color depending on the user
%input
figure
if zColor
    p = patch('Faces',smoothISurface.faces,'Vertices',smoothISurface.vertices,'FaceVertexCData',color,'FaceColor','interp');
    colormap('jet')
    p.EdgeColor = 'none';
    daspect([2 2 1])
    view(3);
    axis tight
    camlight
    lighting gouraud
    title('Z-coloring')
else
    p2 = patch(smoothISurface);
    p2.FaceColor = colorModel;
    p2.EdgeColor = 'none';
    view(3);
    axis tight
    camlight
    lighting gouraud
    title('unicolor');
end
