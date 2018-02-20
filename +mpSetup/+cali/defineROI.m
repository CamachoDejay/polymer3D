function [ ROI ] = defineROI( c_win, ch_c1, ch_c2, im_size )
%DEFINEROI define rois in the cameras for 3D SOFI
%   Detailed explanation goes here

    c_win = round(c_win);
    assert(c_win(1) <= 512,'Problems with channel size');
    c_del = c_win./2;
    ROI = zeros(8,4);
    
    for ch = 1:8
        
        if ch<5
            p1 = round(ch_c1(ch,:) - c_del);            
        else
            p1 = round(ch_c2(ch-4,:) - c_del);            
        end
        
        p2 = p1 + c_win -1;
        
        if p1 < 1
            er = abs(p1) + 1;
            p1 = p1 + er;
            p2 = p2 + er;
        elseif p2 > 2048
            
            er = p2 - im_size(2);
            p2 = p2 - er;
            p1 = p1 - er;
        end
        
        ROI(ch,:) = [p1(1), p1(2), c_win(1), c_win(2)];
        
    end
        
end

