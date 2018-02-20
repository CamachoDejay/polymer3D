function [ Ich, Ich_norm ] = findChInt( data12, in_focus )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = size(data12,4);
I = zeros(8,N);
I(:,:) = mean(mean(data12,2),1);
Inorm = I./ repmat(max(I,[],2),1,N);
bleaching = mean(Inorm,1);
Icor = I./ repmat(bleaching,8,1);
    
Ich = zeros(8,1);
for i = 1:8
    Ich(i) = Icor(i,in_focus(i).frame);
end

Ich_norm = Ich./max(Ich);
    
end

