function [ C, T, Z, IDF, P, F ] = getInfoFromString( s )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

idx1 = strfind(s,'C=');
idx2 = strfind(s,'T=');
idx3 = strfind(s,'Z=');
idx4 = strfind(s,'IFD=');
idx5 = strfind(s,'PlaneCount=');
idx6 = strfind(s,'><UUID');
idx7 = strfind(s,'FileName=');
idx8 = strfind(s,'>urn:uuid');

C   = s(idx1+3:idx2-8);
T   = s(idx2+3:idx3-8);
Z   = s(idx3+3:idx4-3);
IDF = s(idx4+5:idx5-3); 
P   = s(idx5+12:idx6-2);
F   = s(idx7+10:idx8-2);
end

