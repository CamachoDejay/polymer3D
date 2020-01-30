function [ imShifts ] = simpleImShift2( inFocus, cam1, cam2)
%SIMPLEIMSHIFT fast im shift calculation at the pixel resolution
%   spatial cross-correlation algorithm to determine shift of coordinates
    nPlanes = length(inFocus);
    imShifts = zeros(nPlanes-1,2);
    focus = inFocus(1).frame;
    
    imCh1 = cam1(:,:,1,focus);
    % shift must be small thus
    if nPlanes == 2
        maxShift = 300;
    else
        maxShift = 100;
    end%[pixels]
    % generating mask, so I only look at a local max around shift 0 to
    % max_shift
    bw = false(size(imCh1).*2 - [1 1]);
    bw_c = round(size(bw)./2);
    bw(bw_c(1),bw_c(2)) = true;
    se = strel('disk',maxShift,8);
    bw = imdilate(bw,se);
    
    for chIdx = 1:nPlanes-1
        
        idxPlane1 = [inFocus.globalch]==chIdx;
        idxPlane2 = [inFocus.globalch]==chIdx+1;
        
        focus     = round((inFocus(idxPlane1).frame+inFocus(idxPlane2).frame)/2);
        camIdx1 = inFocus(idxPlane1).cam;
        camIdx2 = inFocus(idxPlane2).cam;
        
        if camIdx1==1
            imCh1 = cam1(:,:,inFocus(idxPlane1).ch,focus);
        else
            imCh1 = cam2(:,:,inFocus(idxPlane1).ch,focus);
        end
         
        if camIdx2==1
            imCh2 = cam1(:,:,inFocus(idxPlane2).ch,focus);
        else
            imCh2 = cam2(:,:,inFocus(idxPlane2).ch,focus);
        end
        
        % cross correlation
        res    = normxcorr2(imCh1,imCh2);
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
    % calculate shift for a certain reference plane
    refPlane = round(nPlanes/2);
    tmpImShifts = zeros(nPlanes,2);
    idx = 1:nPlanes;
    idx(idx==4) = [];
    tmpImShifts(idx,:) = imShifts;
    for i = 1: size(tmpImShifts,1)
        
        if i<refPlane
            tmpImShifts(i,:) = -sum(imShifts(i:refPlane-1,:),1);
        elseif i>refPlane
            tmpImShifts(i,:) = +sum(imShifts(refPlane:i-1,:),1);
        else
            %refPlane so we do nothing
        end
        
    end
    %reorder shifts
    imShifts([inFocus.globalch],:) = tmpImShifts;

end

