classdef SRCalibration < handle
    %SRCalibration Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path
        info
        SRCalMovies
        cal2D
        calib
        accuracy
    end
    
    methods
        
        function obj = SRCalibration(path2SRCal,cal2D,info)
            %zCalibration Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = path2SRCal;
            obj.cal2D = cal2D;
            obj.info = info;
        end
        
        function set.path(obj, path)
            assert(ischar(path), 'Path should be given as a string');
            assert(isfolder(path), 'The path given is not a folder, ZCalibration expect a folder. In the folder it is expected to find separate folder for each zCalMovie.')
            
            folderContent = dir(path);
            %Get how many folder are in the main folder
            idx = sum(cellfun(@sum,{folderContent.isdir}));
            %Matlab always store ., .. as folder for relative path so we
            %want to find more than 2 folder in folderContent.
            assert(sum(idx)>2, 'No folder was found in the path given. Expected to find separate folder for each zCalMovie.');
            obj.path = path;
            
        end
        
        function set.cal2D(obj,cal2D)
            assert(isstruct(cal2D),'2D Calibration is expected to be received as a structure');
            assert(isfield(cal2D,'fullPath'),'Missing field "fullPath" in cal2D structure');
            assert(isfield(cal2D,'file'),'Missing field "file" in cal2D structure');
            assert(isfield(cal2D,'info'),'Missing field "info" in cal2D structure');
            
            obj.cal2D = cal2D;
            
        end
        
        function retrieveSRCalMov(obj)
            %we get the zCalibration directory
            folder2Mov = dir(obj.path);
            folder2Mov = folder2Mov(cell2mat({folder2Mov.isdir}));
            %loop through the content of the directory
            for i = 3:size(folder2Mov,1)
                %If element i is a folder
                if folder2Mov(i).isdir
                    %Check if the directory
                    currDir = dir([folder2Mov(i).folder filesep folder2Mov(i).name]);
                    idx = contains({currDir.name}, 'ome.tif');
                    if ~all(idx==0)
                        %we extract z motor position to check if the movie
                        %is indeed a zCalibration (expect zStack)
                        tmp = Core.SRCalMovie([folder2Mov(i).folder filesep folder2Mov(i).name], obj.cal2D,obj.info);
                        tmp.calibrate;
                        [zStep, ~] = tmp.getZPosMotor;
                        %TODO: Check other motor position (we do not want
                        %any other movement here.
                        
                        if zStep ~= 0
                            %if it is we store
                            obj.SRCalMovies.(['SRCal' num2str(i-2)]) = tmp;
                            
                        else
                            %if it is not we throw a warning message as it
                            %might be that many movie are in the main
                            %folder
                            warning(['In ' folder2Mov(i).folder filesep folder2Mov(i).name ' no movement of the Z motor was found, the file is therefore ignored']);
                            
                        end
                    else
                        
                        warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ome.Tif and is therefore ignored']);
                        
                    end
                end
            end
        end
        
        function retrieveSRCalData(obj,detectParam, trackParam)
            nPlanes = obj.SRCalMovies.(['SRCal' num2str(1)]).calibrated.nPlanes;
            %Checking user input
            assert(nargin==3, 'retrieveSRCalData expects 3 inputs, 1)detection Parameters, fit z parameter, tracking parameter');
            assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
            assert(and(isstruct(trackParam),and(isfield(trackParam,'euDistPx'),isfield(trackParam,'commonPlanes'))),'trackParam is expected to be a struct with 2 fields "euDistPx", the tolerated euclidian distance between consecutive frames, "commonPlanes" nPlane to be considered as frame consistent ([1 9])');
            
            %Extraction of SRData
            nfields = numel(fieldnames(obj.SRCalMovies));
            allData = cell(nPlanes-1,1);
            calPerPlane = cell(nPlanes,1);
            for i = 1: nfields
                disp(['Retrieving data from SRCal file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                if i == 1
                    %Ask user for info about the setup for detection
                    obj.SRCalMovies.(['SRCal' num2str(i)]).giveInfo;
                    
                else
                    %get the info about the setup stored into the first
                    %object
                    obj.SRCalMovies.(['SRCal' num2str(i)]).info = obj.SRCalMovies.(['SRCal' num2str(1)]).getInfo;
                    
                end
                %Molecule detection
                obj.SRCalMovies.(['SRCal' num2str(i)]).findCandidatePos(detectParam);
                
                 %SR fitting
                obj.SRCalMovies.(['SRCal' num2str(i)]).SRLocalizeCandidate;
                
                %plane consolidation
                frames = 1:obj.SRCalMovies.(['SRCal' num2str(i)]).calibrated.nFrames;
                obj.SRCalMovies.(['SRCal' num2str(i)]).consolidatePlanes(6,frames,detectParam.consThresh)
                
                %getting Calibration data
                [SRCalibData,dataPerPlane] = obj.SRCalMovies.(['SRCal' num2str(i)]).getSRCalData(trackParam);
                
                %Get data for every particle every planes together stored
                %in allData.
                
                for j = 1:length(allData)
                    allData{j} = [allData{j} ; SRCalibData{j} ];
                    calPerPlane{j} = [calPerPlane{j}; dataPerPlane{j}(dataPerPlane{j}.ellip<1,:)];
                    calPerPlane{j+1} = [calPerPlane{j+1};dataPerPlane{j+1}(dataPerPlane{j+1}.ellip>1,:)];
                end
                
            end
            
            %reorder allData
            
            disp('=================> DONE <===================');
            %store the results of the data extraction
            obj.calib.trackParam = trackParam;
            obj.calib.data = allData;
            obj.calib.dataPerPlane = calPerPlane;
            
        end
        
        function [transMat,corrData] = corrTranslation(obj,refPlane)
            assert(~isempty(obj.calib.data),'You need to extract the SR calibration data before correction for translation');
            SRCalibData = obj.calib.data;
            data2Corr = obj.calib.dataPerPlane;
            %Calculate the translation
            [transMat] = Core.SRCalMovie.getTrans(SRCalibData);
            
            %Correct the data
              %Correct the data for rotation
              corrData = cell(size(data2Corr));
            for i = 1:length(data2Corr)
                currentData = data2Corr{i};
                [cData] = Core.SRCalMovie.applyTrans(currentData,transMat,refPlane);
                corrData{i} = cData;
            end
            
            %storing
            obj.calib.SRCorrData = corrData;
            obj.calib.corr.trans = transMat;
            
            %Saving
            SRCal = obj.calib.corr;
            fileName = sprintf('%s%sSRCalibration.mat',obj.path,'\');
            save(fileName,'SRCal');

        end
        
        function [rotMat,corrData] = corrRotation(obj,refPlane)
            
            SRCalibData = obj.calib.data;
            data2Corr = obj.calib.dataPerPlane;
           %Calculate the rotation
            [rotMat] = Core.SRCalMovie.getRot(SRCalibData);
            
            corrData = cell(length(data2Corr),1);
             %Correct the data for rotation
            for i = 1:length(data2Corr)
                currentData = data2Corr{i};
                [cData] = Core.SRCalMovie.applyRot(currentData,rotMat,refPlane);
                corrData{i} = cData;
            end

            obj.calib.SRCorrData = corrData;
            obj.calib.corr.rot = rotMat;

            %Saving
            SRCal = obj.calib.corr;
            fileName = sprintf('%s%sSRCalibration.mat',obj.path,'\');
            save(fileName,'SRCal');
            
        end
        
        function checkAccuracy(obj,refPlane)
            
            data = obj.calib.dataPerPlane;
            corrData = obj.calib.SRCorrData;
            
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
                %Here we compare the data of plane x with ellipticity <1
                %with the data of plane x+1 with ellipticity >1 before and
                %after correction.
                
                mData = data{i};
                nData = data{i+1};
                
                errRow{i} = mData.row(mData.ellip<1) - nData.row(nData.ellip>1);
                errCol{i} = mData.col(mData.ellip<1) - nData.col(nData.ellip>1);
                euclDist{i} = sqrt(errRow{i}.^2+errCol{i}.^2);
                
                cData = corrData{i};
                dData = corrData{i+1};
                
                errCRow{i}   =  cData.row(cData.ellip<1) - dData.row(dData.ellip>1);
                errCCol{i}   =  cData.col(cData.ellip<1) - dData.col(dData.ellip>1);
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
        
    end
    
    methods (Access = private)
        
        
    end
end

