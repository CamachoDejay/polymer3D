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
        
        function SRCalibrate(obj,trackParam)
            
            %#1 Track particle in Z (= Consolidation between frames)
            [traces, counter] = obj.trackInZ(trackParam);
            
            %#2 Extract Data per particles
            [partData] = obj.extractPartData;
            
            %#3 Find frames where a particle is approx equally defocused in
            %2 consective planes.
            [idx2Frame]= obj.findDefocusedFrame(partData);
            
            %#4 Extract the data per plane
            [SRCalibData] = obj.getCalibData(partData,idx2Frame);
            obj.SRCalData = SRCalibData;
            
            %#5 Find Translation
            [transMat] = obj.getTranslation(SRCalibData);
            obj.SRCorr.translation = transMat;
            
        end
        
        function applyCalib(obj,refPlane)
            disp('Starting super-resolution correction of the data');
            data = obj.SRCalData;
            nPlanes = length(data);
            correction = obj.SRCorr;
            corrData = data;
            for i = 1 : nPlanes
                currentPlane = i;
                if currentPlane < refPlane
                    
                    idx2Corr = currentPlane:refPlane-1;
                    sign = -1;
                    
                elseif currentPlane > refPlane
                     
                    idx2Corr = refPlane:currentPlane-1;
                    sign = +1;
                else
                    idx2Corr = [];
                end
                %1 Translation
                if ~isempty(idx2Corr)
                    for j = 1:length(idx2Corr)
                        corrData{i}(:,1:2) = corrData{i}(:,1:2) + sign* correction.translation{idx2Corr(j)};
                    end
                end
                
            end
            obj.SRCorrData = corrData;
            disp('=============> DONE <==============');
        end
        
        function checkCalib(obj)
           data = obj.SRCalData;
           corrData = obj.SRCorrData;
           err = cell(length(data)-1,1);
           errCorr = err;
           
           for i =1:length(data) - 1
               
               err{i} = data{i,1}(data{i,1}(:,3)<1,1:2) - data{i+1,1}(data{i+1,1}(:,3)>1,1:2);
               errCorr{i} = corrData{i,1}(corrData{i,1}(:,3)<1,1:2) - corrData{i+1,1}(corrData{i+1,1}(:,3)>1,1:2);
           
           end
           
           figure
           fullErr = [];
           fullCErr = [];
           for i =1:length(data) - 1
               
               [errN, errEdge] = histcounts(err{i});
               errBin = errEdge(1:end-1)+(errEdge(2)-errEdge(1));
               [errCN, errCEdge] = histcounts(errCorr{i});
               errCBin = errCEdge(1:end-1)+(errCEdge(2)-errCEdge(1));
               
               
               subplot(1,3,1)
               hold on
               bar(errBin,errN)
               xlim([-2 2])
               title('Before SR Correction');
               xlabel('Error in Pixel');
               ylabel('Occurennce');
               hold off
               subplot(1,3,2)
               hold on
               bar(errCBin,errCN)
               xlim([-2 2])
               title('After SR Correction');
               xlabel('Error in Pixel');
               ylabel('Occurennce');
               hold off
               
               fullErr = [fullErr; err{i}];
               fullCErr = [fullCErr; errCorr{i}];
               
           end
           
               [fErrN, fErrEdge] = histcounts(fullErr);
               fErrBin = fErrEdge(1:end-1)+(fErrEdge(2)-fErrEdge(1));
               [fErrCN, fErrCEdge] = histcounts(fullCErr);
               fErrCBin = fErrCEdge(1:end-1)+(fErrCEdge(2)-fErrCEdge(1));
               
               subplot(1,3,3)
               hold on
               plot(fErrBin,fErrN);
               plot(fErrCBin,fErrCN);
               legend({'Before Corr','After Corr'});
               hold off
           
           avgErr = mean(abs(fullErr));
           stdErr = std(abs(fullErr));
           avgCErr = mean(abs(fullCErr));
           stdCErr = std(abs(fullCErr));
           
           fprintf('The absolute error before correction is: %d \n', mean(avgErr));
           fprintf('The accompanying standard deviation is: %d \n', mean(stdErr));
           fprintf('The absolute error after correction is: %d \n', mean(avgCErr));
           fprintf('The accompanying standard deviation is: %d \n', mean(stdCErr));
           
           
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
        
        function [transMat] = getTranslation(obj,SRCalibData)
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
        
      
        
    end
    
    
end

