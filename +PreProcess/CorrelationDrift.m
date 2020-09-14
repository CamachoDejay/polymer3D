function [correctedStack,Drift]=CorrelationDrift(im_In,scalingFactor,correlationInfo)

    %maxDrift = correlationInfo.maxDrift;
    corrSz = correlationInfo.corrSz; 
    driftPeriod = correlationInfo.driftPeriod;

    frame1 = squeeze(im_In(:,:,1));

    %find center of data
    BW = imbinarize(uint16(frame1));
    cent = regionprops(BW,'centroid');
    tmp = [cent.Centroid];
    coord(:,1) = tmp(1:2:end);
    coord(:,2) = tmp(2:2:end);
    xIdx = round(mean(coord(:,1)));
    yIdx = round(mean(coord(:,2)));

    %Crop the image to save computing time in the correlation
    lb = [xIdx - corrSz + 1, yIdx - corrSz + 1];
    ub = [xIdx + corrSz - 1, yIdx + corrSz - 1];
    %Check that the cropping occurs within the frame
    if lb(1)<1
        lb(1) = 1;
    end

    if lb(2)<1
        lb(2) = 1;
    end

    if ub(1)>size(im_In,1)
        ub(1) = size(im_In,1);
    end

    if ub(2)>size(im_In,2)
        ub(2) = size(im_In,2);
    end

    %cropping
    imCropped = im_In(lb(2): ub(2),lb(1):ub(1),:);
    %scaling for subpixel drift
    if scalingFactor>1
        imCropped = imresize(imCropped,scalingFactor);
        im_In = imresize(im_In,scalingFactor);
    end

    correctedStack = zeros(size(im_In));
    Drift=zeros(size(imCropped,3),2);

    refFrame = squeeze(mean(imCropped(:,:,1:driftPeriod),3));
    for i = 1:floor(size(imCropped,3)/driftPeriod)
        corrMatrix = normxcorr2(refFrame,...
            squeeze(mean(imCropped(:,:,(i-1)*driftPeriod+1:i*driftPeriod),3)));
        [valy,indy] = max(corrMatrix,[],1);
        [~,ind2y]   = max(valy);

        Drift((i-1)*driftPeriod+1:i*driftPeriod,1)  =...
            -ceil(size(corrMatrix,1)/2) + indy(ind2y);

        [valx,indx] = max(corrMatrix,[],2);
        [~,ind2x]   = max(valx);
        Drift((i-1)*driftPeriod+1:i*driftPeriod,2) =...
            -ceil(size(corrMatrix,2)/2) + indx(ind2x);
    end

    for i=1:size(imCropped,3)
        currentIm =im_In(:,:,i);
        [correctedStack(:,:,i)] = PreProcess.correctImageDrift(currentIm,...
            -Drift(i,:));
    end

    %Rescale image
    if scalingFactor>1
        correctedStack = imresize(correctedStack,1/scalingFactor);
        Drift = Drift/scalingFactor;
    end

end