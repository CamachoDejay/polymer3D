function [ im_shifts ] = simpleImShift( in_focus, data12 )
%SIMPLEIMSHIFT fast im shift calculation at the pixel resolution
%   spatial cross-correlation algorithm to determine shift of coordinates
    
    im_shifts = zeros(8,2);
    fr = in_focus(1).frame;
    im_ch1 = data12(:,:,1,fr);
    % shift must be small thus
    max_shift = 100; %[pixels]
    bw = false(size(im_ch1).*2 - [1 1]);
    bw_c = round(size(bw)./2);
    bw(bw_c(1),bw_c(2)) = true;
    se = strel('disk',max_shift,8);
    bw = imdilate(bw,se);
    
    for i = 2:8
    
        fr     = in_focus(i).frame;
        im_chi = data12(:,:,i,fr);
        res    = normxcorr2(im_ch1,im_chi);
        res    = res.*bw;

        [~,m_corr] = max(res(:));
        [mx,my]=ind2sub(size(res),m_corr);
        mx=-mx+size(im_ch1,1);
        my=-my+size(im_ch1,2);
        
        im_shifts(i,:) = [mx my];
        
    end
    
end

