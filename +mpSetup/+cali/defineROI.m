function [ ROI ] = defineROI( chaWin, chC1, chC2, im_size )
%DEFINEROI define ROIs in the cameras for multiplane experiments
    nChan = size(chC1,1);
    chaWin = round(chaWin);
    if (chaWin(1)>= 512)
        warning('Channel size is unexpectedly large if you are using 8 planes configuration');
    end
    chaDel = repmat(chaWin,nChan,1)./2;
    ROI = zeros(nChan*2,4);
    
    cornerC1 = round(chC1 - chaDel);
    cornerC2 = round(chC2 - chaDel);
    ROI(:,1:2) = cat(1,cornerC1,cornerC2);
    ROI(:,3) = chaWin(1);
    ROI(:,4) = chaWin(2);
    x1 = ROI(:,1);
    x2 = x1+chaWin(1)-1;
    y1 = ROI(:,2);
    y2 = y1+chaWin(2)-1;
    assert(and(all(x1>0),all(x2<im_size(2))),'Problems with ROI lims')
    assert(and(all(y1>0),all(y2<im_size(1))),'Problems with ROI lims')
    
%     for ch = 1:8
%         % check which camera im looking at
%         if ch<5
%             % looking at first cam
%             p1 = round(chC1(ch,:) - chaDel);            
%         else
%             %looking at second cam
%             p1 = round(chC2(ch-4,:) - chaDel);            
%         end
%         
%         p2 = p1 + chaWin -1;
%         
%         if p1 < 1
%             er = abs(p1) + 1;
%             p1 = p1 + er;
%             p2 = p2 + er;
%         elseif p2 > 2048
%             
%             er = p2 - im_size(2);
%             p2 = p2 - er;
%             p1 = p1 - er;
%         end
%         
%         ROI(ch,:) = [p1(1), p1(2), chaWin(1), chaWin(2)];
%         
%     end
        
end

