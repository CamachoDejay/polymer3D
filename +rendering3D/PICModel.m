addpath('..\');
%compile c code 
mex +rendering3D\smoothpatch_curvature_double.c -v
mex +rendering3D\smoothpatch_inversedistance_double.c -v
mex +rendering3D\vertex_neighbours_double.c -v
%% load data
path025 = 'N:\Project Z-stacks\New stacks\1K_0.25mg_30 nm_sample 1\tif_1K_0.25 mg-30 nm_Sample 1\SegmentedStacks\Seg_adapt_Sequence 001_Job 2_002_l0_z00.tif';
path05  = 'N:\Project Z-stacks\New stacks\1K_0.5mg_30 nm_sample 1\Tif_1K_0,5mg_30 nm_sample 1\SegmentedStacks\Seg_adapt_Job 6_008.tif';

%%
path2Use = path05;

fInfo = Load.Movie.tif.getinfo(path2Use);
%
myVolume = Load.Movie.tif.getframes(path2Use,1:fInfo.Frame_n);

% cropping
myVol = myVolume(400:565,400:565,:);
myVol = imcomplement(myVol);
%%  

test = isosurface(myVol,1/2);

%% smoothing

test2 = rendering3D.smoothpatch(test,0,100,1);
%%
color = test2.vertices(:,3)/max(test2.vertices(:,3));
uniColor = [0.5137,0.6941,0.9882];

figure
p = patch('Faces',test2.faces,'Vertices',test2.vertices,'FaceVertexCData',color,'FaceColor','interp');
colormap('jet')
%p = patch(test2.vertices(:,1),test2.vertices(:,2),test2.vertices(:,3),uniColor);

%p.FaceColor = uniColor;
p.EdgeColor = 'none';
daspect([2 2 1])
view(3); 
axis tight
camlight 
lighting gouraud