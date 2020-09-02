%% Importing data from .SPE file into matlab array 'a'
% Imports SPE file and returns the stack of frames along with the
% multiplication gain used during acquisition.
%
% Ex:
%
% [I,MG]=ImportSPE;
%


function [A,MG,header]=ImportSPE(folder,file)

    
    
    fid = fopen([folder '/' file],'r'); 

    header = fread(fid,2050,'uint16=>uint16'); % 2050 uint16 = 4100 bytes = 32800 bits

    Xdim = double(header(22));
    Ydim = double(header(329));
    Zdim = double(header(724));
    MG = header(2049);
    DataType = header(55);

    switch DataType
        case 0	% FLOATING POINT (4 bytes / 32 bits)
            ImMat = fread(fid,inf,'float32=>float32');
        case 1	% LONG INTEGER (4 bytes / 32 bits)
            ImMat = fread(fid,inf,'int32=>int32');
        case 2	% INTEGER (2 bytes / 16 bits)
            ImMat = fread(fid,inf,'int16=>int16');
        case 3	% UNSIGNED INTEGER (2 bytes / 16 bits)
            ImMat = fread(fid,inf,'uint16=>uint16');
    end

    fclose(fid);

    a = reshape(ImMat,Xdim,Ydim,Zdim);
    % clear ImMat; %clear some memory

    %permute the X and Y dimensions so that an image looks like in Winview
    a = permute(a,[2,1,3]);
    A = double(a);
    %cd('E:\Chemical Physics')
    [~,file,~]=fileparts(file);
    file=regexprep(file,'\s','');

end