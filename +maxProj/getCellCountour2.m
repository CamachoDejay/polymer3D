function [contour,pContour] = getCellCountour2(BW,IM,thresh1,doplot)
%% get contour of polymer around the cell
BW = ~BW;
%plot raw image
    if doplot
        figure(1)
        subplot(1,4,1)
        imagesc(IM)
        colormap('hot')
        axis image
    end
   
    
    BWL = bwlabel(BW);
    
    %Get the largest area
    cBWarea = regionprops(BW,'Area');
    [val,idx2BiggestArea] = max(cell2mat({cBWarea.Area}));
   % [~,id] = min(val);
    if isempty(idx2BiggestArea)
    else
    %kill all the other area found
        BW(BWL~=idx2BiggestArea) = 0;
    end
    
     % Clean up boundary
    se = strel('disk',10);
    BW = imclose(BW,se);
   
    
    %get outside contour
    [pContour] = bwboundaries(BW);
    
    %find largest contour
    contourSize = cellfun(@size,pContour,'UniformOutput',false);
    [tmp] = cellfun(@max,contourSize,'UniformOutput',false);
    [~,idx2OuterContour] = max(cell2mat(tmp));
    %only keep largest countour
    if isempty(idx2OuterContour)
        pContour = [];
    else
        pContour = pContour{idx2OuterContour};
    end
     if doplot
        figure(1)
        subplot(1,4,2)
        imagesc(IM)
        hold on
       
        plot(pContour(:,2),pContour(:,1),'w');
        colormap('hot')
        axis image
    end
%% get inner countour of the polymer around the cell == contour of the cell

    %make cropbased on contour
    %nIM = IM(min(pContour(:,1)):max(pContour(:,1)),min(pContour(:,2)):max(pContour(:,2)));
    mask  = poly2mask(pContour(:,2),pContour(:,1),size(IM,1),size(IM,2));
    nIM = double(IM).*mask;
    %gBW = imbinarize(nIM,thresh1);
    gBW = nIM;
    gBW(gBW<thresh1*max(max(gBW)))  = 0;
    gBW(gBW>=thresh1*max(max(gBW))) = 1;
    
    se = strel('disk',4);
    gBW = imclose(gBW,se);
    gBW = ~gBW;
    gBW = imclearborder(gBW);
    
    se = strel('disk',8);
    gBW = imclose(gBW,se);    
    
   %%
    %get label
    gBWL = bwlabel(gBW);
    
    %get area
    gBWarea = regionprops(gBW,'Area');
    %keep 2nd biggest area
    [val,idx2BiggestArea] = max(cell2mat({gBWarea.Area}));
    if isempty(idx2BiggestArea)
    else
    %kill all the other area found
        
        gBW(gBWL~=idx2BiggestArea)= 0;
    end

    % Clean up boundary
    se = strel('disk',10);
    gBW = imclose(gBW,se);
      
     if doplot
        figure(1)
        subplot(1,4,3)
        imagesc(gBW)
        axis image
        colormap('hot')
    end
    %get contour
    contour = bwboundaries(gBW);
    %find largest contour
    contourSize = cellfun(@size,contour,'UniformOutput',false);
    [tmp] = cellfun(@max,contourSize,'UniformOutput',false);
    [~,idx2OuterContour] = max(cell2mat(tmp));
    %only keep largest countour
    if isempty(idx2OuterContour)
    contour = [];
    else
        contour = contour{idx2OuterContour}; 
    end
      
     if doplot
        figure(1)
        subplot(1,4,4)
        imagesc(IM)
        hold on
        plot(contour(:,2),contour(:,1),'w','LineWidth',1.5);
        colormap('hot')
        axis image
    end

end