function [contour,pContour] = getCellCountour(IM,thresh1,thresh2,doplot)
%% get contour of polymer around the cell


    if doplot
        figure(1)
        subplot(1,4,1)
        imagesc(IM)
        colormap('hot')
        axis image
    end
        
    %plot raw image
    BW = imbinarize(IM,thresh1);
    BWL = bwlabel(BW);

    %Get the largest area
    cBWarea = regionprops(BW,'Area');
    [~,idx2BiggestArea] = max(cell2mat({cBWarea.Area}));
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
    nIM = IM(min(pContour(:,1)):max(pContour(:,1)),min(pContour(:,2)):max(pContour(:,2)));
    thresh = adaptthresh(nIM,0.8);
%     thresh = multithresh(nIM,4);
%     gBW = imquantize(nIM,thresh);
    %
    gBW = imbinarize(nIM,thresh);
    se = strel('disk',3);
    gBW = imclose(gBW,se);
    gBW = ~gBW;  

 
    
   %%
    %get label
    gBWL = bwlabel(gBW);
    %get area
    gBWarea = regionprops(gBW,'Area');
    %keep 2nd biggest area
    [val,idx2BiggestArea] = maxk(cell2mat({gBWarea.Area}),2);
    if isempty(idx2BiggestArea)
    else
    %kill all the other area found
        [~,idx] = min(val);
        gBW(gBWL~=idx2BiggestArea(idx))= 0;
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
        contour(:,1) = contour(:,1) + min(pContour(:,1)); 
        contour(:,2) = contour(:,2) + min(pContour(:,2)); 
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