


path2File = 'C:\Users\Boris\Desktop\SegmentTest\Stack2Segment\Tif_1K05mg_30nm_sample1\Job 2_004.tif';

frames2Load = 1:152;
IM     = Load.Movie.tif.getframes(path2File, frames2Load); %Loading on of the frame

S = 1;
% size of pixel in z vs x/y
pixZ  = 4;
zFactor = 2;
sigma = [S,S,S*zFactor/pixZ];
IMs = imgaussfilt3(IM, sigma);


path2File = 'C:\Users\Boris\Desktop\SegmentTest\Stack2Segment\Tif_1K05mg_30nm_sample1';
filename = [path2File filesep 'stackGauss.tif'];
dataStorage.nBTiff(filename,IMs,8);