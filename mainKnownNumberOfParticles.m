clear 
clc
close all
%% Create the Movie object

path2File= 'E:\Data\Leuven Data\2018\ZHao\TestCode\400 nm AuNPs 1064 nm laser stepwise - 2 polarization_1';
info.type = 'transmission';
myMov = Core.Movie(path2File,info);
myMov.giveInfo;
%%
myMov.showFrame(1,2)

%% get the data from the movie
fullStack = myMov.getFrame;

%%
for i = 1:nFrames
    
    
    
end

