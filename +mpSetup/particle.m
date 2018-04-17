classdef particle < handle
    %PARTICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frame
        approxPos
        ROI
        image
        superResData
        superResLoc
%         FWHMpx
        gPSFsigPx
    end
    
    methods
        function obj = particle(pos,frame,ROI,data)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.approxPos = pos;
            obj.frame = frame;
            obj.ROI = ROI;
            obj.image = data(ROI(1):ROI(2),ROI(3):ROI(4),:);
        end
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end

        function show(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if length(obj) == 1
                figure()
                for i = 1:8
                    subplot(2,4,i)
                    imagesc(obj.image(:,:,i))
                    axis image
                    title(['Im plane: ' num2str(i)])
                end
            else
                % object is a list
                for objIdx = 1:length(obj)
                    tmpObj = obj(objIdx);
                    tmpObj.show;
                end
                
            end
            
        end
        
        function obj = setPSFprops(obj,objNA, emWave, pxSizeNm)
            % for a water immersion obj the best-fit gasuss to the PSF has
            % sigma = 0.25 wavelength / NA
            if length(obj) == 1
                sigma_nm = 0.25 * emWave/objNA;
                sigma_px = sigma_nm / pxSizeNm;
                obj.gPSFsigPx = sigma_px;
            else
                % object is a list
                for objIdx = 1:length(obj)
                    tmpObj = obj(objIdx);
                    tmpObj.setPSFprops(objNA, emWave, pxSizeNm);
                end
                
            end
                
            
            
        end
        
        function obj = superResolve(obj)
            
            if length(obj) == 1
                data = obj.image;
                sig = [obj.gPSFsigPx, obj.gPSFsigPx];
                assert(~isempty(sig),'cant do fit, you forgot to setPSFprops')
                srOut = nan(8,6);
                % note that due to edge efects on the ROI def the ROI
                % center might not be at obj.approxPos
                ROIcent = [mean(obj.ROI(1:2)), mean(obj.ROI(3:4))];
%                 stdIm = zeros(8,1);
                
                sizeX = size(data,1);
                sizeY = size(data,2);
                % number of pixels in each image
                pixN = sizeX*sizeY;
                % centroid of the window
                xid   = 1:sizeX;
                yid   = 1:sizeY;
                cVal  = [mean(xid), mean(yid)];
                for i = 1:8
                    % load data
                    im = data(:,:,i);
           
                    % get phasor output
                    [x,y,e] = Localization.phasor(im);
                    srOut(i,1) = ROIcent(1)+x;
                    srOut(i,2) = ROIcent(2)+y;
                    srOut(i,3) = e;

                    %   change image to doubles
                    im = double(im);
                    % calculate std
                    srOut(i,4) = std(im(:));
                    
                    % get likelihood ratio test
                    %   get gaussian model
                    %       position
                    pos = [cVal(1)+x,cVal(2)+y];
                    %       function
                    [G] = Misc.gaus2D(pos, sig, xid,yid,1);
                    G = G';
                    
                    %   calculation of G_bar and sGbar2.
                    %       Mean of the gaussian PSF
                    G_mean  = mean(G(:));
                    %       G bar: gaussian PSF minus its mean
                    G_bar = G - G_mean;
                    %       sum of G bar squared
                    sGbar2 = sum(sum(G_bar.^2));
                    %       convolution of the image with G_bar
                    imConv = sum(sum(G_bar.*im));
                    %       calculation of I_hat
                    I_hat  = imConv/sGbar2;
                    %       denominator, which depends on how good it fits
                    %       to random noise
                    den = var(im(:))*pixN;
                    %   calculation of the difference betwen the likelihood
                    %   of having the data explained by random noise or by
                    %   a gaussian. L(H_0)-L(H_1)
                    LRT = ( (pixN)/2 ) * (log(1-((I_hat^2 * sGbar2)/den)));
                    srOut(i,5) = LRT;
                    
                    % get fit to gaussian
                    bg = ones(size(G));
                    % where A are the separate spectrum and b is what what we measured
                    meas   = im(:);
%                     meas   = meas - mean(meas);
%                     meas   = meas ./ std(meas);
                    shapes = [G(:), bg(:)];
                    
                    % find the coefficients that solve the system of equations
                    % [x, resN] = lsqnonneg(shapes,meas);
                    eps = mldivide(shapes,meas);
                    fit = shapes*eps;
                    % fit = reshape(fit,size(im));
                    dif = meas-fit;
                    RMSD = sqrt( sum(dif.^2)/length(dif));
                    srOut(i,6) = RMSD;

                end

                obj.superResData = srOut;
                
                % find imaging plane with best focus using LRT
                [~, bestFocus] = min(srOut(:,5));
                
                % for the moment we store a single number, in the future we
                % will use the mp information
                obj.superResLoc = [srOut(bestFocus,:), bestFocus];
                
                
            else
                % object is a list
                for objIdx = 1:length(obj)
                    tmpObj = obj(objIdx);
                    tmpObj.superResolve();
                end
                
            end
            
            
        end
    end
end

