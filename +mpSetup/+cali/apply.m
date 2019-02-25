function [data] = apply( cam1, cam2, cal )
%APPLY rearranges data and corrects for intensity diffecences
%between channels
%   Detailed explanation goes here
if or(isempty(cam1),isempty(cam2))
    error('Calibration expect 2 cameras, update will be performed later');
else
    h = waitbar(0,'Please wait applying calibration');
    % if need be we flip camera 2, this is generally the case
    if cal.flipCam2
        waitbar(.1,h,'Flipping cam')
        cam2 = flip(cam2,2);
    end
    waitbar(.2,h,'Gettingg channel data')
    [ chC1, chC2 ] = mpSetup.cali.getChData( cam1, cam2, cal.ROI );
    sTmp = size(chC1);
    sTmp(3) = sTmp(3)*2;
    data = ones(sTmp,'uint16');
    
    if cal.correctInt
        C = cal.Icorrf;
    else
        C = uint16(ones(8,1));
        data = uint16(data);
    end
    
    if cal.reorder
        newor = [cal.neworder];
    else
        newor = 1:sTmp(3);
    end
  
    % correct int
    waitbar(.5,h,'Doing some simple math...')
    for i = 1:size(chC1,3)
        data(:,:,i,:) = chC1(:,:,i,:).*C(i);
    end
    
    waitbar(.7,h,'Doing some simple math...')
    for i = 1:size(chC1,3)
        data(:,:,i+size(chC1,3),:) = chC2(:,:,i,:).*C(i+size(chC1,3));
    end
    
    waitbar(.9,h,'Reordering...')
    % reorder planes
    data = data(:,:,newor,:);
    
    close(h)
end
end
