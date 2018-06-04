function [data, frameInfo, movInfo] = loadAndCal( fPath, cal, frames )
%LOADDATA Summary of this function goes here
%   Detailed explanation goes here

    

% load general information about the multi-plane movie
[frameInfo, movInfo, ~ ]= Load.Movie.ome.getInfo( fPath );


switch nargin
    case 2
        disp('Loading all frames lets hope we have enough memory')
        frames = 1:movInfo.maxFrame;
    case 3
        disp('Loading user defined frames')
    otherwise
        error('not enought input arguments')
end

% load the raw data 
[ movC1, movC2] = Load.Movie.ome.load( frameInfo, movInfo, frames );

[data] = mpSetup.cali.apply( movC1, movC2, cal );

end

