function [data] = loadAndCal( fPath, cal )
%LOADDATA Summary of this function goes here
%   Detailed explanation goes here

% load general information about the multi-plane movie
[frameInfo, movInfo, ~ ]= loadMovie.ome.getInfo( fPath );

% load the raw data 
[ movC1, movC2] = loadMovie.ome.load( frameInfo, movInfo, 1:movInfo.maxFrame );

[data] = mpSetup.cali.apply( movC1, movC2, cal );

end

