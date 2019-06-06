function [ C, T, Z, IDF, P, F, pos, expT ] = getInfoFromString( tifStr, planeStr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%TiffData
idx1 = strfind(tifStr,'C=');
idx2 = strfind(tifStr,'T=');
idx3 = strfind(tifStr,'Z=');
idx4 = strfind(tifStr,'IFD=');
idx5 = strfind(tifStr,'PlaneCount=');
idx6 = strfind(tifStr,'><UUID');
idx7 = strfind(tifStr,'FileName=');
idx8 = strfind(tifStr,'>urn:uuid');

%PlaneData
idx9  = strfind(planeStr,'TheC=');
idx10 = strfind(planeStr,'TheT=');
idx11 = strfind(planeStr,'TheZ=');
idx12 = strfind(planeStr,'PositionX=');
idx13 = strfind(planeStr,'PositionXUnit=');
idx14 = strfind(planeStr,'PositionY=');
idx15 = strfind(planeStr,'PositionYUnit=');
idx16 = strfind(planeStr,'PositionZ=');
idx17 = strfind(planeStr,'PositionZUnit=');
idx18 = strfind(planeStr,'ExposureTime=');
idx19 = strfind(planeStr,'ExposureTimeUnit=');

C   = tifStr(idx1+3:idx2-8);
C2  = planeStr(idx9+6:idx10-3);
assert(strcmp(C,C2),'Unexpected C')

T   = tifStr(idx2+3:idx3-8);
T2  = planeStr(idx10+6:idx11-3);
assert(strcmp(T,T2),'Unexpected T')

Z   = tifStr(idx3+3:idx4-3);
Z2  = planeStr(idx11+6:end-1);
assert(strcmp(Z,Z2),'Unexpected Z')

IDF = tifStr(idx4+5:idx5-3); 
P   = tifStr(idx5+12:idx6-2);
F   = tifStr(idx7+10:idx8-2);

Xpos   = str2double(planeStr(idx12+11:idx13-3));
Ypos   = str2double(planeStr(idx14+11:idx15-3));
Zpos   = str2double(planeStr(idx16+11:idx17-3));

pos = [Xpos, Ypos, Zpos];

expT   = str2double(planeStr(idx18+14:idx19-3));
end

