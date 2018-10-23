function [ focus_met, in_focus, fit ] = getFocusMetric( chData1c, chData2c, Z1, Z2 )
%GETFOCUSMETRIC gets information about when the channels are in focus.
%Inspired in the work done in EPFL for 3D SOFI

    N = size(chData1c,4);
    nChan = size(chData1c,3);
    focus_met = zeros(N, nChan*2);
    fit = zeros(N, nChan*2);
    in_focus(nChan*2).cam   =  [];
    in_focus(nChan*2).ch    =  [];
    in_focus(nChan*2).frame =  [];
    in_focus(nChan*2).zpos  =  [];
    
    for i = 1:nChan
        in_focus(i).cam   =  1;
        in_focus(i).ch    =  i;
        % get image sequence for 2 channels to compare
        tmp=squeeze(chData1c(:,:,i,:));
        
        % see when each channel is in focus
        focus_met(:,i) = squeeze(mean(max(tmp)));
        [~,in_focus(i).frame] = max(focus_met(:,i));
        [ zFocus, fit(:,i) ] = mpSetup.cali.getSubResPlanePosition(focus_met(:,i),Z1);
        in_focus(i).zpos = zFocus;
    end
    
    for i = nChan+1:nChan*2
        in_focus(i).cam   =  2;
        in_focus(i).ch    =  i-nChan;
        % get image sequence for 2 channels to compare
        tmp=squeeze(chData2c(:,:,i-nChan,:));
        
        % see when each channel is in focus
        focus_met(:,i) = squeeze(mean(max(tmp)));
        [~,in_focus(i).frame] = max(focus_met(:,i));
        [ zFocus, fit(:,i) ] = mpSetup.cali.getSubResPlanePosition(focus_met(:,i),Z2);
        in_focus(i).zpos = zFocus;
    end
    
end

