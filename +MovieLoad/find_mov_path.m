function [ mov_list ] = find_mov_path( init_dir, known_extensions )
%FIND_MOV_PATH Summary of this function goes here
%   Detailed explanation goes here

% path to tiff file, regardless of it is a stack or not  
[main_folder] = uigetdir(init_dir, 'Select the directory that contains all movies to be analyzed');

if ~main_folder
    disp('No path given, boombing out...')
    mov_list = [];
    return
end

mov_path = dir([main_folder filesep '*.*']); 
    
% list of all files
ct = sum(~[mov_path.isdir]);
mov_list = cell(ct,1);

c=1;
for Index = 1:length(mov_path)
    if ~mov_path(Index).isdir
        baseFileName = mov_path(Index).name;
        baseDir = [main_folder filesep baseFileName];
        
        [~, ~, extension] = fileparts(baseDir);
                
        isKnown = sum( strcmp( known_extensions,extension)) == 1;

        if isKnown
            mov_list{c} = [main_folder filesep baseFileName];
            c=c+1;
        else
            warning((['I dont know/trust this file extension: ' extension]))
        end
    end
end

mov_list(c:end)=[];
    
end

