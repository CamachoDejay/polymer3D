function [ focus_met, in_focus ] = getFocusMetric( chData1c, chData2c, Z1, Z2 )
%GETFOCUSMETRIC gets information about when the channels are in focus.
%Inspired in the work done in EPFL for 3D SOFI

    N = size(chData1c,4);
    focus_met = zeros(N, 8);
    
    in_focus(8).cam   =  [];
    in_focus(8).ch    =  [];
    in_focus(8).frame =  [];
    in_focus(8).zpos  =  [];
    
    for i = 1:4
        in_focus(i).cam   =  1;
        in_focus(i).ch    =  i;
        % get image sequence for 2 channels to compare
        tmp=squeeze(chData1c(:,:,i,:));
        
        % see when each channel is in focus
        focus_met(:,i) = squeeze(mean(max(tmp)));
        [~,in_focus(i).frame] = max(focus_met(:,i));
        in_focus(i).zpos = Z1(in_focus(i).frame);
    end
    
    for i = 5:8
        in_focus(i).cam   =  2;
        in_focus(i).ch    =  i-4;
        % get image sequence for 2 channels to compare
        tmp=squeeze(chData2c(:,:,i-4,:));
        
        % see when each channel is in focus
        focus_met(:,i) = squeeze(mean(max(tmp)));
        [~,in_focus(i).frame] = max(focus_met(:,i));
        in_focus(i).zpos = Z2(in_focus(i).frame);
    end
    
end

