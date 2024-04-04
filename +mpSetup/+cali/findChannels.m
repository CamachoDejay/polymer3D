function [ chC, bgC, common_w ] = findChannels( im, doFigure,nChan )
%FINDCHANNELS finds where the channels are in the calibration files for
%multiplane setup. We take as input an time -average or -max image of a
%camera and then use simple integration to find the areas of fluorescence
%and background
    %From Rafa:
    sig = im(20:end-20,:);
    sig = sig(:);
    bg = im(:,1000:1040);
    bg = bg(:);
    
    tHold = Misc.tholdSigBg(bg,sig);
    im = im>tHold;
    %remove small pixel
    im = bwareaopen(im,21);
    % Clean up boundary
    se = strel('square',5);
    im = imclose(im,se);
    im = bwareaopen(im,16);
    
    switch nargin
        case 1
            doFigure = false;
    end
    
    assert(size(im,2)==2048, 'unexpected movie size, are you working with the sCMOS?')
    % integration over the cols, here we expect to find 4 channels
    s1 = sum(im,1);
    % find change points between fluorescence and bg
    chP = findCp(s1,'bottom',nChan);
    % now we want to find the size of the window and its center, also that
    % of the bg
    x1 = chP(:,1);
    x2 = chP(:,2);
    % channels x width
    chXw = x2-x1 +1;
    % channels center point
    chC = zeros(size(chP,1),2);
    chC(:,1) = x1 + chXw./2;
    % width and center of the bg
    bgW  = x1(2:end)-x2(1:end-1);
    bgC  = round(x2(1:end-1) + (bgW./2));
    % minimun window size that can be used to set the ROI.
    commonXw = min(chXw);

    % now we look at the different channels in the rows. We know that there
    % are small vertical shift that we whish to compensate in a rowgh way.
    chLims = [1; bgC; size(im,2)];
    chYw = zeros(size(chP,1),1);
    for i = 1:size(chP,1)
        chIm = im(:,chLims(i):chLims(i+1));
        % integration
        s2 =  sum(chIm,2);
        % change points
        chP2 = findCp(s2,'bottom',1);
        % windown size
        chYw(i) =  chP2(2) - chP2(1) + 1;
        chC(i,2) = chP2(1) + chYw(i)/2;
    end
    
    
    % min win size that can be used to set the ROI
    commonYw = min(chYw);
    % storing the window sized used
    common_w = [commonXw commonYw];
    
    if doFigure
        corner = round(chC - [chXw, chYw]./2);
        figure()
        imagesc(im)
        axis image
        for i = 1:size(chP,1)
            rectangle('Position', [corner(i,1) corner(i,2) chXw(i), chYw(i)])
        end
    end
end

function ch_p = findCp(trace_in,sCase,nCP)
    % nCP number of expected change points
    trace_in = trace_in(:);
    % smoothing the trace
    sT = smooth(trace_in,20);
    % normalizing the trace
    sT = sT - min(sT);
    sT = sT ./ max(sT);
    upCP = nan;
    doCP = nan;
    % clearing edge effects
    sT(1:10) = min(sT);
    sT(end-9:end) = min(sT);
    % choosing the staring point
    switch sCase
        case 'bottom'
            thVal = 0.09;
            dTh   = 0.01;
            
        case 'top'
            thVal = mean(sT);
            dTh   = -0.01;
            
        otherwise
            error('dont know what to dom8, either top or bottom')
    end
    full = false;
    go = true;
    % iterative thresholding up to the point that we get the right results
    while go
        % threshold value
        thVal = thVal+dTh;
        if or(thVal<=0, thVal>=1)
            warning('could not find the channels');
            full = true;
            break;
        end
        
        % binarized trace
        binT = sT > (thVal);
        if strcmp(sCase,'bottom')
%             disp('ok')
            binT(1) = 0;
            binT(end) = 0;
        end
        % forcing the ends to be false - this is important so we alway have
        % at least a CP at the firts and last row
%         binT(1) = false;
%         binT(end) = false;
        dYl   = diff(binT);
        upCP  = find(dYl==1);
        doCP  = find(dYl==-1);
        if and(~isempty(upCP), ~isempty(doCP))
            tmp1 = upCP(1);
            tmp2 = doCP(end);
            doCP(doCP < tmp1) = [];
            upCP(upCP > tmp2) = [];
        end
        nUp   = length(upCP);
        nDo   = length(doCP);
        
        
        if nUp==nDo
            wSize = doCP-upCP;
            if and(nUp==nCP, all(wSize>0))
                % we are finish, so we stop
                if all(wSize>200)
                    go = false;
                end
                
            end
        end

    end
    % we delete data up to 10 so if the change point happen at 10 most
    % likely there was no change point.
    if  upCP == 10
        upCP =0;
    end
    if doCP == length(trace_in) - 10
        doCP = length(trace_in)-1;
    end
    
    ch_p = [upCP+1, doCP+1];
    if full
        ch_p = [1 length(trace_in)];
    end
end

