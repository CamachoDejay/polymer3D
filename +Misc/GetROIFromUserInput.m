function [Mask] = GetROIFromUserInput(MPMovie,ROI)
%     imshow(double(MPMovie(:,:,1))./max(double(MPMovie(:,:,1))))
    nChan = 1:size(ROI,1);
    movC1 = MPMovie;
    Mask = zeros(size(movC1,1),size(movC1,2));

    for i=1:nChan(length(nChan)/2)
            col1 = ROI(i,1); 
            col2 = ROI(i,1) + ROI(i,3) - 1;
            row1 = ROI(i,2);
            row2 = ROI(i,2) + ROI(i,4) - 1;
            newimg = max(MPMovie(row1:row2,col1:col2,:),[],3);
            if i==1
                figure()

                imshow(mat2gray(newimg));
                title('Please, select region of interest by drawing a rectangle around it')

                h = imrect;
                Pos = h.getPosition;
                x0 =Pos(1);
                y0 =Pos(2);
                len = Pos(3);
                wid = Pos(4);
            end
            newimg(:,:)=0;
            newimg(floor(y0):floor(y0)+floor(wid) , floor(x0):floor(x0)+floor(len),:)=1;
            
            Mask(row1:row2,col1:col2,1) = newimg;
    end
end