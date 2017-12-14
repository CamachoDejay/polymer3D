function [ em_pos ] = getRandPos( size, emitter_number )
%GETRANDPOS random positions of emitters with the contrain that the image
%is square
%   Detailed explanation goes here

em_n = emitter_number;
em_pos = rand(em_n,2).*size;
end

