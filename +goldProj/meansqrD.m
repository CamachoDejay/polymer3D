function [D MSD]=meansqrD(cod)

%[D MSD]=meansqrD(cod, stp)
% the function gives the distance ('D') between each point and
% mean-squre-dispracement ('MSD')
%
%input parameters; cod is 2D aray including coordinate of trajectry.


%%%Fist calculate distance between each frame
D=[];
idx=1:length(cod);
DX=diff(cod(idx,1).'); DY=diff(cod(idx,2).');
D=[D sqrt(DX.^2 + DY.^2)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cnt=1;MSD=[];
%Calculate mean-squere-displacement
for i=1:length(cod)-1
    stp=i;cnt=1;D1=[];D2=[];DX=[];DY=[];;
    while cnt<=stp & cnt+stp<=length(cod)
        idx=cnt:stp:length(cod);
        DX=diff(cod(idx,1).');
        DY=diff(cod(idx,2).');
        D1=[D1 sqrt(DX.^2 + DY.^2)];
        cnt=cnt+1;
        if ~isempty(D1)
            D2=D1(~isnan(D1));
            if ~isempty(D2)
                MSD(i)=mean(D2.^2);
            else
                MSD(i)=NaN;
            end
        end
    end %while
end

D=D';
MSD=MSD';