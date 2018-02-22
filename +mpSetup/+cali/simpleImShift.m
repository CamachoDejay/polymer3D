function [ imShifts ] = simpleImShift( inFocus, cam1, cam2)
%SIMPLEIMSHIFT fast im shift calculation at the pixel resolution
%   spatial cross-correlation algorithm to determine shift of coordinates
    
    imShifts = zeros(8,2);
    focus = inFocus(1).frame;
    
    imCh1 = cam1(:,:,1,focus);
    % shift must be small thus
    max_shift = 100; %[pixels]
    % generating mask, so I only look at a local max around shift 0 to
    % max_shift
    bw = false(size(imCh1).*2 - [1 1]);
    bw_c = round(size(bw)./2);
    bw(bw_c(1),bw_c(2)) = true;
    se = strel('disk',max_shift,8);
    bw = imdilate(bw,se);
    
    for chIdx = 2:8
    
        focus     = inFocus(chIdx).frame;
        if chIdx<5
            % I look at first cam
            imChi = cam1(:,:,chIdx,focus);
        else
            % I look at second cam
            imChi = cam2(:,:,chIdx-4,focus);
        end
        % cross correlation
        res    = normxcorr2(imCh1,imChi);
        % applying the mask
        res    = res.*bw;
        % finding shift
        [~,m_corr] = max(res(:));
        [mx,my]=ind2sub(size(res),m_corr);
        mx=-mx+size(imCh1,1);
        my=-my+size(imCh1,2);
        % storing the shifts
        imShifts(chIdx,:) = [mx my];
        
    end
    
end

