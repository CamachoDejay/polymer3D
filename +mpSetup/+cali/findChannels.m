function [ ch_c, bg_c, common_w ] = findChannels( ave_im )
%FINDCHANNELS finds where the channels are in the calibration files for
%multiplane setup

    assert(size(ave_im,2)==2048, 'unexpected movie size')
    s1 = sum(ave_im,1);
    ch_p = find_cp(s1,'bottom',4);
    x1 = ch_p(:,1);
    x2 = ch_p(:,2);
    
    ch_xw = x2-x1 +1;
    ch_c  = x1 + ch_xw./2;
    
    bg_w  = x1(2:end)-x2(1:end-1);
    bg_c  = round(x2(1:end-1) + (bg_w./2));
    
%     t_w = [bg_c; 2048] - [1; bg_c];
    
    common_xw = min(ch_xw);
%%
    ch1 = ave_im(:,1:bg_c(1));
    s2 =  sum(ch1,2);
    ch_p = find_cp(s2,'top',1);
    ch1_y1 = ch_p(:,1);
    ch1_y2 = ch_p(:,2);
    ch_yw(1) =  ch1_y2 - ch1_y1 + 1;
    ch_c(1,2) = ch1_y1 + ch_yw(1)/2;
    
    ch2 = ave_im(:,bg_c(1):bg_c(2));
    s2 =  sum(ch2,2);
    ch_p = find_cp(s2,'top',1);
    ch2_y1 = ch_p(:,1);
    ch2_y2 = ch_p(:,2);
    ch_yw(2) =  ch2_y2 - ch2_y1 + 1;
    ch_c(2,2) = ch2_y1 + ch_yw(2)/2;
    
    ch3 = ave_im(:,bg_c(2):bg_c(3));
    s2 =  sum(ch3,2);
    ch_p = find_cp(s2,'top',1);
    ch3_y1 = ch_p(:,1);
    ch3_y2 = ch_p(:,2);
    ch_yw(3) =  ch3_y2 - ch3_y1 + 1;
    ch_c(3,2) = ch3_y1 + ch_yw(3)/2;
    
    ch4 = ave_im(:,bg_c(3):end);
    s2 =  sum(ch4,2);
    ch_p = find_cp(s2,'top',1);
    ch4_y1 = ch_p(:,1);
    ch4_y2 = ch_p(:,2);
    ch_yw(4) =  ch4_y2 - ch4_y1 + 1;
    ch_c(4,2) = ch4_y1 + ch_yw(4)/2;
    
    common_yw = min(ch_yw);
    
    common_w = [common_xw common_yw];
    ch_c = round(ch_c);
    

end

function ch_p = find_cp(y_in,sCase,nCP)
    % nCP number of expected change points
    y_in = y_in(:);
    yy = smooth(y_in,20);
    yy = yy - min(yy);
    yy = yy ./ max(yy);
    upCP = nan;
    doCP = nan;
    switch sCase
        case 'bottom'
            thVal = 0.0;
            dTh   = 0.01;
            
        case 'top'
            thVal = mean(yy);
            dTh   = -0.01;
            
        otherwise
            error('dont know what to dom8, either top or bottom')
    end
    
    go = true;
    while go
        thVal = thVal+dTh;
        yl = yy > (thVal);
        yl(1) = false;
        yl(end) = false;
        dYl = diff(yl);
        upCP = find(dYl==1);
        doCP = find(dYl==-1);
        go = ~and(length(upCP)==nCP, length(doCP)==nCP);
        if or(thVal<=0, thVal>=1)
            error('could not find the channels')
        end
    end
    
    ch_p = [upCP+1, doCP+1];
end

