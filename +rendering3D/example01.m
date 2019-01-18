% this is an example from mathworks
clear
close all
clc

% load example of volumetric data from MATLAB:
%   D: contains MRI data of a human head - uint8
%   map: gray map 89x3 double
%   siz: size of the MRI data 128,128,27. Meaning 27 slides of 128x128 pix

load mri
% we use squeeze to remove the singleton dimension the data
D = squeeze(D);

%% vis the data as 2D slides
figure(1)
for i = 1:size(D,3)
    imagesc(D(:,:,i))
    title(sprintf('Slice %d',i))
    axis image
    drawnow
    pause(0.2)
end

%% vis 3D contour slices
figure(2)
cm = brighten(jet(length(map)),-.5);
colormap(cm)
slice2see = [1,12,19,27];
contourslice(D,[],[],slice2see,8);
view(3);
axis tight

%% to open volumeViewer use
volumeViewer(D)

%%
clear
V=flow;
figure(1)
for i = 1:size(V,3)
    imagesc(V(:,:,i))
    title(sprintf('Slice %d',i))
    axis image
    drawnow
    pause(0.2)
end

%%
figure(2)
slice2see = [6,12,18];
contourslice(V,[],[],slice2see,8);
view(3);
axis tight
%% 
% generating 3D volumetric data
[x_grid,y_grid,z_grid,vals] = flow;
% we want to generate a surface from the data. For this we will use
% isosurface, which conects pixels in the volume of the same value.
%   value to look for in the data and conect
val2conect = -3;
%   generation of the isosurface
isoS = isosurface(x_grid, y_grid, z_grid, vals, val2conect);
%   creation of filled poligon for display
pS   = patch(isoS);
%   smooth the rendering via isonormals
isonormals(x_grid,y_grid,z_grid,vals,pS)
%   modifing some properties of the poligon
%       color of the faces
pS.FaceColor = 'red';
%       transparency of the faces
pS.FaceAlpha = 0.5;
%       color of the edges - mesh, in this case we remove it
pS.EdgeColor = 'none';

% hold becuase I want to create and plot another isofurface
hold on

% see comments above for details
val2conect = 0;
isoS2 = isosurface(x_grid, y_grid, z_grid, vals, val2conect);
pS2   = patch(isoS2);
isonormals(x_grid,y_grid,z_grid,vals,pS2)
pS2.FaceColor = 'green';
pS2.FaceAlpha = 0.5;
pS2.EdgeColor = 'none';

% data aspect ratio for the current axes. The data aspect ratio is the
% relative length of the data units along the x-axis, y-axis, and z-axis.
daspect([1 1 1])
% setting the view angle - 3 is the default three-dimensional view, az =
% ?37.5, el = 30
view(3)
% axis limits equal to the data range
axis tight
% create a camera light object, this is just to make the rendering look
% more realistic
camlight
% now we specify how to generate the lights gouraud calculates the vertex
% normals and interpolates linearly across the faces
lighting gouraud