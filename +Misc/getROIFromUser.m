function [Mask] = getROIFromUser(im)
    
    Mask = zeros(size(im,1),size(im,2)); 
    %Show image to user
    figure()
    imagesc(im);
    colormap('gray');
    axis image
    title('Please, select region of interest by drawing a rectangle around it')
    %allow user to select ROI
    h   = imrect;
    Pos = wait(h);
    x0  = floor(Pos(1));
    y0  = floor(Pos(2));
    len = floor(Pos(3));
    wid = floor(Pos(4));
   
    Mask(y0:y0+wid-1,x0:x0+len-1) = 1;

end