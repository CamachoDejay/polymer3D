classdef SRCalMovie < Core.ZCalMovie
    %SuperResCal will hold information about the operation that needs to be
    %performed on the coordinates to match them to the plane of reference.
    
    %Although it is not obvious why it should inherit from zCalMovie, the
    %reason is that in essence, the superResCalMovie is a zStack and
    %similar analysis (e.g. localization, consolidation, trackInZ,...) will
    %need to be performed prior to calculate the superResCal.
    
    %Thus rather than having two classes with a lot in common (thus
    %duplicate) it was decided that it would inherit from zCalMovie
    
    properties (SetAccess = 'private')
        SRCalData
        SRCorr
        SRCorrData
    end
    
    methods
        %Constructor
        function obj = SRCalMovie(raw, cal)
            
            obj  = obj@Core.ZCalMovie(raw,cal);
            
        end
        
        function [SRCalibData] = getSRCalData(obj,trackParam)
            switch nargin
                case 1
                    error('Need tracking parameters to SR-Calibrate')
                case 2
                otherwise
                    error('Unexpected number of input');
            end
            disp('Starting SR data extraction')
            %#1 Track particle in Z (= Consolidation between frames)
            [~,~] = obj.trackInZ(trackParam);
            
            %#2 Extract Data per particles
            [partData] = obj.extractPartData;
            
            %#3 Find frames where a particle is approx equally defocused in
            %2 consective planes.
            [idx2Frame]= obj.findDefocusedFrame(partData);
            
            %#4 Extract the data per plane
            [SRCalibData] = obj.getCalibData(partData,idx2Frame);
            obj.SRCalData = SRCalibData;
            disp('==========> DONE ! <============');
        end
        
        function [transMat,corrData] = corrTranslation(obj,refPlane)
            assert(~isempty(obj.SRCalData),'You need to extract the SR calibration data before correction for translation');
            SRCalibData = obj.SRCalData;
            %Calculate the translation 
            [transMat] = obj.getTrans(SRCalibData);
            
            %Correct the data
            [corrData] = obj.corrTrans(SRCalibData,transMat,refPlane);
            
            obj.SRCorrData = corrData;
        end
        
        function [transMat,rotMat,corrData] = corrRotation(obj,refPlane)
            
            [transMat,corrData] = obj.corrTranslation(refPlane);
            
             %Calculate the rotation
            [rotMat] = obj.getRot(corrData);
            
            %Correct the data
            [corrData] = obj.corrRot(corrData,rotMat, refPlane);
   
            obj.SRCorrData = corrData;
            
        end
        
        function checkAccuracy(obj,refPlane)
            
            data = obj.SRCalData;
            corrData = obj.SRCorrData;
            
            errRow = cell(length(data)-1,1);
            errCol = errRow;
            euclDist = errRow;
            errCRow = errRow;
            errCCol = errRow;
            euclCDist= errRow;
           
            allDist = [];
            allCorrDist = [];
            %Calculate the error on euclidian distance before and after
            %correction for each individual beads used for calibration
            for i = 1:length(data)-1
                
                    errRow{i} = data{i,1}(data{i,1}(:,3)<1,1) - data{i+1,1}(data{i+1,1}(:,3)>1,1);
                    errCol{i} = data{i,1}(data{i,1}(:,3)<1,2) - data{i+1,1}(data{i+1,1}(:,3)>1,2);
                    euclDist{i} = sqrt(errRow{i}.^2+errCol{i}.^2);
                    
                    errCRow{i} = corrData{i,1}(corrData{i,1}(:,3)<1,1) - corrData{i+1,1}(corrData{i+1,1}(:,3)>1,1);
                    errCCol{i} = corrData{i,1}(corrData{i,1}(:,3)<1,2) - corrData{i+1,1}(corrData{i+1,1}(:,3)>1,2);
                    euclCDist{i} = sqrt(errCRow{i}.^2+errCCol{i}.^2);

                    allDist = [allDist; data{i,1}(:,1) data{i,1}(:,2)];                   
                    allCorrDist = [allCorrDist; corrData{i,1}(:,1) corrData{i,1}(:,2)];
                    
            end
            
            %Plot histogram of errors
            figure
            fullErrRow = [];
            fullCErrRow = [];
            errorPerPlane = zeros(length(data)-1,2);
            errorCorrPerPlane = zeros(length(data)-1,2);
            for i = 1:length(data)-1
                [errN, errEdge] = histcounts(errRow{i},[-2:0.1:2]);
                errBin = errEdge(1:end-1)+(errEdge(2)-errEdge(1));
                [errCN, errCEdge] = histcounts(errCRow{i},[-2:0.1:2]);
                errCBin = errCEdge(1:end-1)+(errCEdge(2)-errCEdge(1));
                
                subplot(1,3,1)
                hold on
                bar(errBin,errN)
                xlim([-2 2])
                title('Before SR Correction');
                xlabel('X Error in Pixel');
                ylabel('Occurennce');
                legend;
                hold off
                subplot(1,3,2)
                hold on
                bar(errCBin,errCN)
                xlim([-2 2])
                title('After SR Correction');
                xlabel('X Error in Pixel');
                ylabel('Occurennce');
                legend;
                hold off
               
                fullErrRow = [fullErrRow; errRow{i}];
                fullCErrRow = [fullCErrRow; errCRow{i}];
                
                errorPerPlane(i,:) = [mean(abs(errRow{i})) mean(errRow{i})];
                errorCorrPerPlane(i,:) = [mean(abs(errCRow{i})) mean(errCRow{i})];
                
            end
            
            [fErrN, fErrEdge] = histcounts(fullErrRow,[-2:0.1:2]);
            fErrBin = fErrEdge(1:end-1)+(fErrEdge(2)-fErrEdge(1));
            [fErrCN, fErrCEdge] = histcounts(fullCErrRow,[-2:0.2:2]);
            fErrCBin = fErrCEdge(1:end-1)+(fErrCEdge(2)-fErrCEdge(1));
            %Compare all errors before and after
            subplot(1,3,3)
            hold on
            plot(fErrBin,fErrN);
            plot(fErrCBin,fErrCN);
            legend({'Before Corr','After Corr'});
            title('Distribution before/after')
            xlim([-2 2])
            hold off
            
            figure
            plane = 1:7;
            height = min(errorPerPlane(:,2)):0.01: max(errorPerPlane(:,2));
            x = ones(1,length(height))*refPlane;
            subplot(1,2,1)
            hold on
            plot(plane,errorPerPlane(:,1));
            plot(plane,errorCorrPerPlane(:,1));
            plot(x,height,'k--')
            legend({'Before Correction','After Correction', 'Ref plane'});
            xlabel('Planes')
            ylabel('error compare to reference plane')
            hold off ;
            
            subplot(1,2,2)
            hold on
            plot(plane,errorPerPlane(:,2));
            plot(plane,errorCorrPerPlane(:,2));
            plot(x,height,'k--')
            legend({'Before Correction','After Correction', 'Ref plane'});
            xlabel('Planes')
            ylabel('error compare to reference plane')
            hold off ;
            
            %display some meaningful values to know what is happening
            avgAbsErr = mean(abs(nonzeros(fullErrRow)));
            avgErr = mean(nonzeros(fullErrRow));
            stdErr = std(abs(nonzeros(fullErrRow)));
            avgAbsCErr = mean(abs(nonzeros(fullCErrRow)));
            avgCErr = mean(nonzeros(fullCErrRow));
            stdCErr = std(abs(nonzeros(fullCErrRow)));
            
            fprintf('The absolute error before correction is: %0.3f pixel\n', mean(avgAbsErr));
            fprintf('The error before correction is: %d \n', mean(avgErr));
            fprintf('The accompanying standard deviation is: %0.3f pixel \n', mean(stdErr));
            fprintf('The absolute error after correction is: %0.3f pixel\n', mean(avgAbsCErr));
            fprintf('The error after correction is: %d \n', mean(avgCErr));
            fprintf('The accompanying standard deviation is: %0.3f pixel\n', mean(stdCErr));
            
            %Let us check that the correction worked or display warning
            if or(avgAbsErr<avgAbsCErr, stdErr<stdCErr)
                warning(['Mean error or standard deviation was smaller before'...
                    ' correction, something might have gone wrong...'])
            end
            
            figure
            subplot(1,3,1)
            plot(allDist(:,2),allDist(:,1),'k+')
            subplot(1,3,2)
            plot(allCorrDist(:,2),allCorrDist(:,1),'r+');
            subplot(1,3,3)
            hold on
            plot(allCorrDist(:,2),allCorrDist(:,1),'r+');
            plot(allDist(:,2),allDist(:,1),'k+')
            hold off 
        end
    end
    
    methods (Access = private)
        
        function [partData] = extractPartData(obj)
            list = obj.particles.List;
            traces = obj.traces.trace;
            nTraces = obj.traces.nTrace;
            
            %Allocate Memory
            partData = cell(1,nTraces);
            %Extraction of particule data occurs here
            for i = 1 : length(traces)
                if ~isempty(traces{i})
                    for j = 1 : length(traces{i})
                        idx2Part = traces{i}{j};
                        if ~isnan(idx2Part)
                            startIdx = size(partData{idx2Part},1)+1;
                            currentData = list{i}{j}(~isnan(list{i}{j}(:,1)),[1:3, end]);
                            idx = repmat(i,size(currentData,1),1);
                            partData{idx2Part}(startIdx:startIdx+size(currentData,1)-1,:) = [list{i}{j}(~isnan(list{i}{j}(:,1)),[1:3, end]) idx];
                        else
                        end
                    end
                end
            end
            
            % Here we clean the data keeping only the particles that were
            % observed in every planes
            for i = 1: size(partData,2)
                
                data2Test = partData{i};
                test = unique(data2Test(:,4));
                
                if length(test)<obj.calibrated.nPlanes
                    %put empty cells where test fails
                    partData{i} = [];
                end
            end
            %Deleting empty cells of the cell array
            partData(cellfun('isempty',partData)) = [];
            
        end
        
        function [defocusFrame] = findDefocusedFrame(obj,partData)
            
            defocusFrame = cell(size(partData));
            
            for i = 1:size(partData,2)
                currentData = partData{i};
                planes = unique(currentData(:,4));
                for j = 1:max(planes)-1
                    dataPlaneA = currentData(currentData(:,4) == planes(j),:);
                    dataPlaneB = currentData(currentData(:,4) == planes(j+1),:);
                    
                    idx = obj.findOptimalDefocusing(dataPlaneA,dataPlaneB);
                    
                    defocusFrame{i}(j,:) = [j idx];
                    
                end
            end
        end
        
        function [idx] = findOptimalDefocusing(obj,dataPlaneA,dataPlaneB)
            nFrames = unique(dataPlaneA(dataPlaneA(:,3)<1,end));
            bestVal = [2 1];
            for i = 1: length(nFrames)
                currentFrame = nFrames(i);
                currentVal = [dataPlaneA(dataPlaneA(:,5)==currentFrame,3) dataPlaneB(dataPlaneB(:,5)==currentFrame,3)];
                cmpVal = [abs(1-1/currentVal(1)) abs(1-currentVal(2))];
                
                if  and(abs(cmpVal(1) - cmpVal(2))< abs(bestVal(1) - bestVal(2)), currentVal(2)>1)
                    bestVal = cmpVal;
                    idx = nFrames(i);
                end
                
            end
            
        end
        
        function [SRCalibData] = getCalibData(obj,partData,idx2Frame)
            nPlanes = obj.calibrated.nPlanes;
            SRCalibData = cell(nPlanes,1);
            for i =1:size(partData,2)
                currentData = partData{i};
                for j = 1:nPlanes-1
                    %Data Plane x
                    idx = and(currentData(:,4)==j,currentData(:,5)==idx2Frame{i}(j,2));
                    SRCalibData{j} = [SRCalibData{j}; currentData(idx,:) ];
                    %Data Plane x+1
                    idx = and(currentData(:,4)==j+1,currentData(:,5)==idx2Frame{i}(j,2));
                    SRCalibData{j+1} = [SRCalibData{j+1}; currentData(idx,:)];
                end
            end
            
        end
                
        function [transMat] = getTrans(obj,SRCalibData)
            nPlanes = obj.calibrated.nPlanes;
            transMat = cell(nPlanes-1,1);
            for i = 1:nPlanes-1
                idx2FrameA = SRCalibData{i}(:,3) < 1;
                idx2FrameB = SRCalibData{i+1}(:,3) > 1;
                
                %Center of mass col row for plane a (ellip>1)
                colCMa = mean(SRCalibData{i}(idx2FrameA,1));
                rowCMa = mean(SRCalibData{i}(idx2FrameA,2));
                
                %Center of mass col row plane b (ellip<1)
                colCMb = mean(SRCalibData{i+1}(idx2FrameB,1));
                rowCMb = mean(SRCalibData{i+1}(idx2FrameB,2));
                
                %if we want to add the correction (not subtract)
                colTrans = colCMa-colCMb;
                rowTrans = rowCMa-rowCMb;
                
                %store
                transMat{i} = [colTrans rowTrans];
                
            end
            
            
        end
        
        function [rotMat]   = getRot(obj,SRCalibData)
            nPlanes = obj.calibrated.nPlanes;
            rotMat = cell(nPlanes-1,1);
            
            for i = 1:nPlanes-1
                idx2FrameA = SRCalibData{i}(:,3) < 1;
                idx2FrameB = SRCalibData{i+1}(:,3) > 1;
                
                planeA = SRCalibData{i}(idx2FrameA,1:3);
                planeB = SRCalibData{i+1}(idx2FrameB,1:3);
                
                planeA(:,3) = 0;
                planeB(:,3) = 0;
                
                planeA = planeA - mean(planeA,1);
                planeB = planeB - mean(planeB,1);
                %check that format is okay
                if and(size(planeA,1)~=3, size(planeA,2)==3)
                    planeA = planeA';
                    planeB = planeB';
                elseif size(planeA,1)==3
                else
                    error('Something is wrong with the dimension of your data, expect a three 3D vector (x;y;z)')
                end
                
                %Calculate rotation matrix (Procedure found at
                %http://nghiaho.com/?page_id=671)
                
                %Calculate covariance matrix
                H = planeA*planeB';
                
                %Single value decomposition
                [U, S, V] = svd(H);
                
                %getting the rotation
                rotMat{i} = V*U';
                
%                 %test if okay
%                 euclDist = mean(sqrt((planeA(1,:) - planeB(1,:)).^2 +...
%                     (planeA(2,:) - planeB(2,:)).^2 + ...
%                     (planeA(2,:) - planeB(2,:)).^2));
%                 
%                 rotPlaneA = rotMat{i}*planeA;
%                 euclDistCorr = mean(sqrt((rotPlaneA(1,:) - planeB(1,:)).^2 +...
%                     (rotPlaneA(2,:) - planeB(2,:)).^2 + ...
%                     (rotPlaneA(3,:) - planeB(3,:)).^2));
%                 
%                 Test = euclDistCorr < euclDist;
            end
        end
        
        function [corrData] = corrTrans(obj,corrData,corr,refPlane)%Refplane
            
            nPlanes = length(corrData);
            %Loop through the planes
            for i = 1 : nPlanes
                 currentPlane = i;
                %act depending on whether the current plane is smaller or
                %bigger than the user-selected reference plane
                if currentPlane < refPlane
                    
                    idx2Corr = currentPlane:refPlane-1;
                    sign = -1;
                    
                elseif currentPlane > refPlane
                    
                    idx2Corr = refPlane:currentPlane-1;
                    idx2Corr = fliplr(idx2Corr);
                    sign = +1;
                    
                else
                    
                    idx2Corr = [];
                    
                end
                %1 Translation
                if ~isempty(idx2Corr)
                    for j = 1:length(idx2Corr)
                        corrData{i}(:,1:2) = corrData{i}(:,1:2) + sign* corr{idx2Corr(j)};
                    end
                end
%                     %remove CM 
%                     corrData{i}(corrData{i}(:,3)<1,1:2) = corrData{i}(corrData{i}(:,3)<1,1:2)-mean(corrData{i}(corrData{i}(:,3)<1,1:2));
%                     corrData{i}(corrData{i}(:,3)>1,1:2) = corrData{i}(corrData{i}(:,3)>1,1:2)- mean(corrData{i}(corrData{i}(:,3)>1,1:2));
%                      
            end
            
        end
        
        function [corrData] = corrRot(obj,corrData,corr,refPlane)
            
            nPlanes = length(corrData);
            %Loop through the planes
            for i = 1 : nPlanes
                currentPlane = i;
                %act depending on whether the current plane is smaller or
                %bigger than the user-selected reference plane
                if currentPlane < refPlane
                    
                    idx2Corr = currentPlane:refPlane-1;
                    sign = false;
                    
                elseif currentPlane > refPlane
                    
                    idx2Corr = refPlane:currentPlane-1;
                    idx2Corr = fliplr(idx2Corr);
                    sign = true;
                    
                else
                    
                    idx2Corr = [];
                    
                end
                
                %Rotation
                if ~isempty(idx2Corr)
                    
                    data2Corr = corrData{i}(:,1:3);
                    data2Corr(:,3) = 0;
                    
                    %corrDA = data2Corr(corrData{i}(:,3)<1,:);
                    %corrDB = data2Corr(corrData{i}(:,3)>1,:);
                    CM = mean(data2Corr);
                   % CMb = mean(corrDB);
                        
%                     corrDA = corrDA - CMa;
%                     corrDB = corrDB - CMb;
                    data2Corr = data2Corr - CM;
                    
                    if and(size( data2Corr,1)~=3, size( data2Corr,2)==3)
                         data2Corr =  data2Corr';
                        %corrDB = corrDB';
                    elseif size(planeA,1)==3
                    else
                        error('Something is wrong with the dimension of your data, expect a three 3D vector (x;y;z)')
                    end
                    
                    for j = 1:length(idx2Corr)
                        
                        if sign
                            rot = corr{idx2Corr(j)}';
                        else
                            rot = corr{idx2Corr(j)};
                        end

%                         if isempty(corrDB)
                            
                            data2Corr = (rot*data2Corr);
                            data2Store = data2Corr';
                            corrData{i}(:,1:2) = data2Store(:,1:2)+CM(1:2);
                            
% %                         elseif isempty(corrDA)
%                             
%                             corrDB = (rot*corrDB);
%                             data2Store = corrDB';
%                             corrData{i}(corrData{i}(:,3)>1,1:2) = data2Store(:,1:2)+CMb(1:2);
%                             
%                         else
%                             
%                             corrDA = (rot*corrDA);
%                             corrDB = (rot*corrDB);
%                             data2Store = corrDA';
%                             corrData{i}(corrData{i}(:,3)<1,1:2) = data2Store(:,1:2)+CMa(1:2);
%                             data2Store = corrDB';
%                             corrData{i}(corrData{i}(:,3)>1,1:2) = data2Store(:,1:2)+CMb(1:2);
                            
                        %end
                    end
                end
            end
            
        end
        
    end
end

