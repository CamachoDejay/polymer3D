function [ mov ] = getframes( path2file, frames )
%GETFRAMES Summary of this function goes here
%   Detailed explanation goes here

assert(min(size(frames))==1,'frames must be a vector of positive integers')

f_n = length(frames);

tObj = Tiff(path2file,'r');
l    = tObj.getTag(256);
w    = tObj.getTag(257);
tObj.setDirectory(frames(1));

im1  = tObj.read;
nClass = class(im1);
mov = zeros(w,l,f_n,nClass);



for i = 1:f_n
    f_i = frames(i);
    tObj.setDirectory(f_i)
    movTmp = tObj.read;  
    mov(:,:,i) = movTmp(:,:,1);    
end
tObj.close


end

