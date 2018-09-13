classdef SRCalMovie < Core.MPCalMovie
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
        calDataPerPlane
        
    end
    
    methods
        %Constructor
        function obj = SRCalMovie(raw, cal)
            
            obj  = obj@Core.MPCalMovie(raw,cal);
            
        end
        
        function [SRCalibData, dataPerPlane] = getSRCalData(obj,trackParam)
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
            [SRCalibData, dataPerPlane] = obj.getCalibData(partData,idx2Frame);
            obj.SRCalData = SRCalibData;
            obj.calDataPerPlane = dataPerPlane;
            disp('==========> DONE ! <============');
        end
        
        function [transMat,corrData] = corrTranslation(obj,refPlane)
            assert(~isempty(obj.SRCalData),'You need to extract the SR calibration data before correction for translation');
            SRCalibData = obj.SRCalData;
            data2Corr = obj.calDataPerPlane;
   
            %Calculate the translation 
            [transMat] = obj.getTrans(SRCalibData);
            
            %Correct the data
            corrData = cell(length(data2Corr),1);
             %Correct the data for rotation
            for i = 1:length(data2Corr)
                currentData = data2Corr{i};
                [cData] = obj.applyTrans(currentData,transMat,refPlane);
                corrData{i} = cData;
            end
                      
            obj.SRCorrData = corrData;
            obj.SRCorr.trans = transMat;
        end
        
        function [rotMat,corrData] = corrRotation(obj,refPlane)
            data2Corr = obj.calDataPerPlane;
            
            SRCalibData = obj.SRCalData;
             %Calculate the rotation
            [rotMat] = obj.getRot(SRCalibData);
            
            corrData = cell(length(data2Corr),1);
             %Correct the data for rotation
            for i = 1:length(data2Corr)
                currentData = data2Corr{i};
                [cData] = obj.applyRot(currentData,rotMat,refPlane);
                corrData{i} = cData;
            end

            obj.SRCorrData = corrData;
            obj.SRCorr.rot = rotMat;
            
        end
        
        function checkAccuracy(obj,refPlane)
            
            data = obj.calDataPerPlane;
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
                    mData = data{i,1};
                    nData = data{i+1,1};
                    
                    errRow{i} = mData.row(mData.ellip<1) - nData.row(nData.ellip>1);
                    errCol{i} = mData.col(mData.ellip<1) - nData.col(nData.ellip>1);
                    euclDist{i} = sqrt(errRow{i}.^2+errCol{i}.^2);
                    
                    cData = corrData{i,1};
                    dData = corrData{i+1,1};
                    errCRow{i} =  cData.row(cData.ellip<1) - dData.row(dData.ellip>1);
                    errCCol{i} =  cData.col(cData.ellip<1) - dData.row(dData.ellip>1);
                    euclCDist{i} = sqrt(errCRow{i}.^2+errCCol{i}.^2);

                    allDist = [allDist; mData.row mData.col];                   
                    allCorrDist = [allCorrDist; cData.row cData.col];
                    
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
        
        function correctionUnitTesting(obj,angle,shift,refPlane)
 %/////////////WARNING This function is made to test that the correction of
              %translation and rotation are working fine. It actually modify
              %the data in the object so only use it to test and when you
              %know what you are doing !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            assert(and(angle<360,angle>0),'Angle is expected to be between 0 and 360')
            assert(and(max(size(shift))==3,min(size(shift))==1),'Shift is expected to be a 3D vector')
            
            SRCalDat = obj.SRCalData;
            dataPerPlane = obj.calDataPerPlane;
            
            angleRad = angle*pi/180;
            %Let us modify the Calibration data so the shift and the
            %rotation is known and correspond to the user input
            
            %calculate the rotation matrix
            rotMat = [cos(angle), -sin(angle), 0 ;
                      sin(angle),  cos(angle), 0 ;
                          0            0       1];
            %take care of the Data
            cData = dataPerPlane{1}(:,{'row','col'});
            cData.z = zeros(size(cData.row,1),1);
            cData = table2array(cData);
            for i = 1:length(dataPerPlane)-1
                cDataPrev = cData;
                
                cCM = mean(cData);%Center of mass
                
                %remove the CM and transpose
                cData = (cData - cCM)';
                
                %do the rotation
                cData = rotMat*cData;
                
                %do the translation and transpose back
                cData = (cData+cCM'-shift')';
                
                dataPerPlane{i}(dataPerPlane{i}.ellip<1,{'row','col'}) = array2table(cDataPrev(:,1:2));
                dataPerPlane{i+1}(dataPerPlane{i+1}.ellip>1,{'row','col'}) =array2table(cData(:,1:2));
                
                SRCalDat{i}(SRCalDat{i}.ellip<1,1:2) =  array2table(cDataPrev(:,1:2));
                SRCalDat{i}(SRCalDat{i}.ellip>1,1:2) =  array2table(cData(:,1:2));
            end
                  
            obj.SRCalData = SRCalDat;
            obj.calDataPerPlane = dataPerPlane;
            
            obj.corrTranslation(refPlane);
            obj.checkAccuracy(refPlane);
            
            obj.corrRotation(refPlane);
            obj.checkAccuracy(refPlane);
            
        end
        
        
    end
    
    methods (Static)
        
        function [corrMat]  = getTrans(SRCalibData)
            nL = length(SRCalibData);
            
            transfo = {'1->2';'2->3';'3->4';'4->5';'5->6';'6->7';'7->8'};
            corrMat = table(zeros(nL,1),zeros(nL,1),zeros(nL,1),zeros(nL,1),...
                zeros(nL,1),zeros(nL,1),zeros(nL,1),transfo,'VariableNames',...
                {'rowTrans','colTrans','rowCMa','colCMa','rowCMb','colCMb',...
                'rot','transformation'});
            
            for i = 1:nL
                currentData = SRCalibData{i};
                %split in two plane
                idx2FrameA = currentData.plane == i;
                idx2FrameB = currentData.plane == i+1;
                
                %Center of mass col row for plane a (ellip>1)
                colCMa = mean(currentData.col(idx2FrameA));
                rowCMa = mean(currentData.row(idx2FrameA));
                
                %Center of mass col row plane b (ellip<1)
                colCMb = mean(currentData.col(idx2FrameB));
                rowCMb = mean(currentData.row(idx2FrameB));
                
                %if we want to add the correction (not subtract)
                colTrans = colCMa-colCMb;
                rowTrans = rowCMa-rowCMb;
                
                %store
                corrMat(i,{'rowTrans','colTrans'}) = table(rowTrans, colTrans);
                corrMat(i,{'rowCMa','colCMa'}) = table(rowCMa,colCMa);
                corrMat(i,{'rowCMb','colCMb'}) = table(rowCMb,colCMb);
            end
            
            
        end
        
        function [rotMat]   = getRot(SRCalibData)
            nL = length(SRCalibData);
            rotMat = table(cell(nL,1),cell(nL,1),cell(nL,1),'VariableNames',{'CMa','CMb','rot'});
            for i = 1:nL
                planeA =[];
                planeB =[];
                
                currCalData = SRCalibData{i};
                %split in two planes
                idx2FrameA = currCalData.plane == i;
                idx2FrameB = currCalData.plane == i+1;
                
                planeA(:,1) = currCalData.row(idx2FrameA);
                planeA(:,2) = currCalData.col(idx2FrameA);

                planeB(:,1) = currCalData.row(idx2FrameB);
                planeB(:,2) = currCalData.col(idx2FrameB);
                
                 %Center of mass col row for plane a (ellip>1)
                rowCMa = mean(planeA(:,1));
                colCMa = mean(planeA(:,2));
                CMa = [rowCMa colCMa];
                
                %Center of mass col row plane b (ellip<1)
                rowCMb = mean(planeB(:,1));
                colCMb = mean(planeB(:,2));
                
                CMb = [rowCMb colCMb];
                               
                planeA = planeA - CMa;
                planeB = planeB - CMb;
               
                planeA = planeA';
                planeB = planeB';
               
                
                %Calculate rotation matrix (Procedure found at
                %http://nghiaho.com/?page_id=671)
                
                %Calculate covariance matrix
                H = planeA*planeB';
                
                %Single value decomposition
                [U, S, V] = svd(H);
                
                %getting the rotation
                rotMat(i,{'CMa'}) = table({CMa});
                rotMat(i,{'CMb'})= table({CMb});
                rotMat(i,{'rot'}) = table({V*U'});
                
           end
        end
        
        function [corrData] = applyTrans(data,corr,refPlane)
            %Refplane
                currentPlane = data.plane(1);
                corrData = data;
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
                row = corrData.row;
                col = corrData.col;
                if ~isempty(idx2Corr)
                    for j = 1:length(idx2Corr)
                        
                        row = row +sign*corr.rowTrans(idx2Corr(j));
                        col = col +sign*corr.colTrans(idx2Corr(j));

                    end
                    corrData(:,'row') = table(row);
                    corrData(:,'col') = table(col);
                    
                end
 
            
        end
        
        function [corrData] = applyRot(data,corr,refPlane)
                
                data2Corr = [];
                currentPlane = data.plane(1);
                corrData = data;
                %act depending on whether the current plane is smaller or
                %bigger than the user-selected reference plane
                if currentPlane < refPlane
                    
                    idx2Corr = currentPlane:refPlane-1;
                    sign = true;
                    
                elseif currentPlane > refPlane
                    
                    idx2Corr = refPlane:currentPlane-1;
                    idx2Corr = fliplr(idx2Corr);
                    sign = false;
                    
                else
                    
                    idx2Corr = [];
                    
                end
                
                %Rotation
                if ~isempty(idx2Corr)
                    
                    data2Corr(:,1) = data.row;
                    data2Corr(:,2) = data.col;
                                       
                    for j = 1:length(idx2Corr)
                        
                        if sign
                            rot = corr.rot{idx2Corr(j)};
                            CM2Sub = corr.CMa{idx2Corr(j)};
                            CM2add = corr.CMb{idx2Corr(j)};
                            
                            
                        else
                            rot = corr.rot{idx2Corr(j)}';
                            
                            CM2Sub = corr.CMb{idx2Corr(j)};
                            CM2add = corr.CMa{idx2Corr(j)};
                            

                        end

                            data2Corr = data2Corr - CM2Sub;
                            data2Corr = data2Corr';
                            data2Corr = (rot*data2Corr);
                            data2Corr = data2Corr' + CM2add;

                    end
                    
                    corrData.row = data2Corr(:,1);
                    corrData.col = data2Corr(:,2);
                    
                else
                    corrData = data;
                    
                end
                
                
            
            
        end
        
    end
    
    methods (Access = private)
        
        function [partData] = extractPartData(obj)
            list = obj.particles.List;
            traces = obj.particles.traces;
            nTraces = obj.particles.nTraces;
            
            %Allocate Memory
            partData = cell(1,nTraces);
            %Extraction of particule data occurs here
            for i = 1 : length(traces)
                if ~isempty(traces{i})
                    for j = 1 : length(traces{i})
                        idx2Part = traces{i}{j};
                        if ~isnan(idx2Part)
                            startIdx = size(partData{idx2Part},1)+1;
                            currentData = list{i}{j}(~isnan(list{i}{j}.row),{'row','col','ellip','plane'});
                            idx = table(repmat(i,size(currentData.row,1),1),'VariableNames',{'frame'});
                            partData{idx2Part}(startIdx:startIdx+size(currentData.row,1)-1,:) = [currentData idx];
                        else
                        end
                    end
                end
            end
            
            % Here we clean the data keeping only the particles that were
            % observed in every planes
            for i = 1: size(partData,2)
                
                data2Test = partData{i};
                test = unique(data2Test.plane);
                
                if length(test)<obj.calibrated.nPlanes
                    %put empty cells where test fails
                    partData{i} = [];
                else%we check if data above and below focus is encountered
                    del = false;
                    for j = 1: length(test)
                        data2Test2 = data2Test(data2Test.plane==j,:);
                        if or(isempty(data2Test2(data2Test2.ellip<1,:)),isempty(data2Test2(data2Test2.ellip>1,:)))
                            del = true;
                        end
                    end
                    if del
                        partData{i} = [];
                    end
                end
                
               
            end
            %Deleting empty cells of the cell array
            partData(cellfun('isempty',partData)) = [];
            
        end
        
        function [defocusFrame] = findDefocusedFrame(obj,partData)
            
            defocusFrame = cell(size(partData));
            
            for i = 1:size(partData,2)
                currentData = partData{i};
                planes = unique(currentData.plane);
                for j = 1:max(planes)-1
                    dataPlaneA = currentData(currentData.plane == planes(j),:);
                    dataPlaneB = currentData(currentData.plane == planes(j+1),:);
                    
                    idx = obj.findOptimalDefocusing(dataPlaneA,dataPlaneB);
                    
                    defocusFrame{i}(j,:) = [j idx];
                    
                end
            end
        end
        
        function [idx] = findOptimalDefocusing(obj,dataPlaneA,dataPlaneB)
            nFrames = unique(dataPlaneA.frame(dataPlaneA.ellip<1));
            bestVal = [2 1];
            for i = 1: length(nFrames)
                currentFrame = nFrames(i);
                currentVal = [dataPlaneA.ellip(dataPlaneA.frame==currentFrame) dataPlaneB.ellip(dataPlaneB.frame==currentFrame)];
                cmpVal = [abs(1-1/currentVal(1)) abs(1-currentVal(2))];
                
                if  and(abs(cmpVal(1) - cmpVal(2))< abs(bestVal(1) - bestVal(2)), currentVal(2)>1)
                    bestVal = cmpVal;
                    idx = nFrames(i);
                end
                
            end
            
        end
        
        function [SRCalibData,dataPerPlane] = getCalibData(obj,partData,idx2Frame)
            nPlanes = obj.calibrated.nPlanes;
            SRCalibData = cell(nPlanes-1,1);
            dataPerPlane = cell(nPlanes,1);
            
            for i =1:size(partData,2)
                currentData = partData{i};
                for j = 1:nPlanes-1
                    %Data Plane x
                    idx = and(currentData.plane==j,currentData.frame==idx2Frame{i}(j,2));
                    idx2LEllip = currentData.ellip<=1;
                    idx2HEllip = currentData.ellip>1;
                    idx2LData  = logical(idx.*idx2LEllip);
                    SRCalibData{j} = [SRCalibData{j}; currentData(idx2LData,:) ];
                    %Data Plane x+1
                    idx = and(currentData.plane==j+1,currentData.frame==idx2Frame{i}(j,2));
                    idx2HData = logical(idx.*idx2HEllip);
                    SRCalibData{j} = [SRCalibData{j}; currentData(idx2HData,:)];
                    
                    %data per plane
                    dataPerPlane{j} = [dataPerPlane{j}; currentData(idx2LData,:)];
                    dataPerPlane{j+1} = [dataPerPlane{j+1}; currentData(idx2HData,:)];
                end
            end
            
        end
                        
       
        
        
        
    end
end

