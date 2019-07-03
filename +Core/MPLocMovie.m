classdef MPLocMovie < Core.MPParticleMovie
    %MPLocMovie hold all the method and info necessary for localization,
    %expension of this object will be first centered around Tracking but
    %it will also be easy to extend to STORM/PALM.
    
    properties (SetAccess = 'protected')
        
        SRCal
        ZCal
        corrected
    end
    
    methods
        
        function obj = MPLocMovie(raw, MPCal,info, SRCal, zCal)
            
            obj  = obj@Core.MPParticleMovie(raw,MPCal,info);
            
            switch nargin
                
                
                case 4
                    obj.SRCal
                    obj.ZCal = [];
                case 5
                    
                    obj.SRCal = SRCal;
                    obj.ZCal = zCal;
               
                otherwise
                    error('Too many input arguments');
            end
        end
        
        function set.SRCal(obj,SRCal)
            
            if ~isempty(SRCal)
                assert(isfolder(SRCal), 'The given path is not a folder');

                %Check Given path
                [file2Analyze] = Core.Movie.getFileInPath(SRCal,'SRCalibration.mat');

                if isempty(file2Analyze)
                    error('No SR calibration file found in the given folder');
                else
                    fileName = [file2Analyze.folder filesep file2Analyze.name];
                    cal = load(fileName);
                    field = fieldnames(cal);
                    cal = cal.(field{1});
                    assert(and(isstruct(cal), and(isfield(cal,'trans'),isfield(cal,'rot'))),...
                        'SR calibration is supposed to be a struct with 2 fields');

                    obj.SRCal = cal; 
                end
            else
                obj.SRCal = [];
            end
        end
        
        function set.ZCal(obj,zCal)
            if ~isempty(zCal)
                assert(isfolder(zCal), 'The given path is not a folder');

                %Check Given path
                [file2Analyze] = Core.Movie.getFileInPath(zCal,'zCalibration.mat');

                if isempty(file2Analyze)
                    error('No z calibration file found in the given folder');
                else
                    fileName = [file2Analyze.folder filesep file2Analyze.name];
                    cal = load(fileName);
                    field = fieldnames(cal);
                    cal = cal.(field{1});
                    assert(isstruct(cal),'zCalibration is supposed to be in cells format');
                    assert(and(isfield(cal,'fitZParam'),isfield(cal,'calib')),...
                        'Something is wrong in the fields of your Z calibration');

                    obj.ZCal = cal; 
                end
            else
                obj.ZCal = [];
            end
        end
        
        function applyCorr(obj,rot,refPlane)
            
            %apply SRCal
            obj.applySRCal(rot,refPlane);
            
            %transform ellipticity into Z
            obj.applyZCal;
            
        end
        
        function applySRCal(obj, rot, refPlane)
            assert(~isempty(obj.unCorrLocPos),'You need to find candidate and SR Localized them before applying corrections');
           
            if isempty(obj.corrLocPos)
                
                    obj.corrLocPos = obj.unCorrLocPos;
                    
            end
            
            if(~isempty(obj.SRCal))
            
                if nargin <2
                    refPlane = 5;
                end

                data = obj.unCorrLocPos;
                
                nPlanesCal = size(obj.SRCal.trans,1)+1;
                nPlanesFile = obj.calibrated.nPlanes;
                assert(nPlanesCal == nPlanesFile,'Mismatch between number of planes in SR calibration and file');
            
                
                disp(['Applying SR calibration...']);
                for i = 1 : length(data)
                    currData = data{i};
                    if ~isempty(currData)
                        currPlanes = unique(currData.plane);
                        for j = 1 : length(currPlanes)
                            currentPlane = currPlanes(j);
                            data2Corr = currData(currData.plane==currentPlane,{'row','col','plane'});

                            if rot
                                corrMat = obj.SRCal.rot;
                                [corrData] = Core.MPSRCalMovie.applyRot(data2Corr, corrMat,refPlane);

                            else
                                corrMat = obj.SRCal.trans;
                                 [corrData] = Core.MPSRCalMovie.applyTrans(data2Corr,corrMat,refPlane);                    
                            end

                            %we store the corrected data
                            obj.corrLocPos{i}(currData.plane==currentPlane,{'row','col','plane'}) = corrData;

                        end
                    end

                end
                obj.corrected.XY = true;
                disp('========> DONE ! <=========');
            else
                obj.corrected.XY = false;
                disp('========> DONE ! <=========');
                warning('SR Calibration not found, no correction was applied');
            end

        end
        
        function applyZCal(obj)
            disp('Applying Z Calibration... ');
            assert(~isempty(obj.unCorrLocPos),'Need to fit before applying the calibration');
            if isempty(obj.ZCal)
                
                warning('Z Calibration needed to correct the data, using Intensity instead');
                obj.addInfo('zMethod','Intensity');
                obj.corrected.Z = false;
                disp('========> DONE ! <=========');
            end
            
            if isempty(obj.corrLocPos)
                obj.corrLocPos = obj.unCorrLocPos;
                warning('Z calibration is currently being applied on non-SRCorrected (X-Y) data');
            end
            
            data = obj.corrLocPos; 
            zCal = obj.ZCal;
            zMethod = obj.info.zMethod;
            
            if strcmp(zMethod,'Intensity')
                obj.corrected.Z = false;
               
            elseif strcmp(zMethod,'PSFE')
            
                %we check which method is best:
                [method] = obj.pickZFitMethod;
                
                %Here we translate ellipticity into z position based on
                %calibration
                nPlanesCal = size(zCal.calib,1);
                nPlanesFile = obj.calibrated.nPlanes;
                assert(nPlanesCal == nPlanesFile,'Mismatch between number of planes in Z calibration and file');

                disp('Applying Z Calibration using PSFE and ZCal');
                for i = 1 : length(data)
                    currData = data{i};
                    nPos = size(currData,1);

                    for j = 1 : nPos

                        currentEllip = currData.ellip(j);
                        currentPlane = currData.plane(j);
                        [zPos] = obj.getZPosition(currentEllip,zCal,currentPlane,method);

                        obj.corrLocPos{i}.z(j) = zPos;
                    end
                    
                end
                         %Here we translate the ellipticity range into zRange for each
                %plane

                ellipRange = zCal.fitZParam.ellipRange;
                nPlanes = obj.calibrated.nPlanes;
                zRange = cell(nPlanes,1);
                
                for i = 1 : nPlanes
                    zRange{i} = obj.getZRange(ellipRange,zCal,i,method);
                end
                
                obj.corrected.Z = true;
                obj.calibrated.zRange = zRange;
                
            else
                error('Unknown Z method');
            end
            
            disp('=======> DONE ! <========');
        end
        
        function [locPos] = getLocPos(obj,frames)
             %Extract the position of the candidate of a given frame
            [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
            locPos = obj.corrLocPos{idx};
            
            if isempty(locPos)
                
                warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                
            end
        end
                
        function superResolve(obj)
            disp('super resolving positions ... ');
            data2Resolve = obj.particles.List;
            
            SRList = [];
            for i = 1:length(data2Resolve)
            
                frameData = data2Resolve{i};
                frameData2Store = table(zeros(size(frameData)),...
                    zeros(size(frameData)),zeros(size(frameData)),zeros(size(frameData)),...
                    zeros(size(frameData)),zeros(size(frameData)),zeros(size(frameData)),...
                    zeros(size(frameData)),zeros(size(frameData)),'VariableNames',...
                    {'row','col','z','rowM','colM','zM','intensity','SNR','t'});
                if strcmpi(obj.info.zMethod,'Intensity')
                    fData = obj.getFrame(i);
                end
                for j = 1:length(frameData)
                   
                    partData = frameData{j};
                    
                    switch obj.info.zMethod
                        case 'Intensity'
                            [data] = obj.resolveXYZInt(partData(:,{'row','col','z','ellip','plane'}),fData);
                        case 'PSFE'
                            [data] = obj.resolveXYZ(partData(:,{'row','col','z','ellip','plane'}));
                    end
                        
                    frameData2Store(j,{'row','col','z','rowM','colM','zM'}) = data;
                    frameData2Store.intensity(j) = partData.intensity(3);
                    frameData2Store.SNR(j) = partData.SNR(3);
                    frameData2Store.t(j) = i;
                
                end
             SRList = [SRList;frameData2Store];   
                
            end
            
            %clean up the list
            SRList(isnan(SRList.row),:) = [];
                
            obj.particles.SRList = SRList;    
            disp('========> DONE ! <=========');
            
        end
                   
        function showCorrLoc(obj,frames)
             part = obj.particles.SRList;
            switch nargin
                case 1
                    frames = min(part.t):max(part.t);
                case 2 
                    [frames] = obj.checkFrame(frames,obj.raw.maxFrame(1));
            end
           
           
            figure()
            hold on
            
            sizeMarker =5;
            scatter3(part.col,part.row,part.z,sizeMarker,part.z,'filled')
            axis ij;

            title('all Localization plotted');
            xlabel('x position in nm');
            ylabel('y position in nm');
            zlabel('z position in nm');
            
            
            hold off
            
        end
        
    end
    
    methods (Static)
        
                    
    end
    
    
    methods (Access = protected)
        
        function [corrData] = applyTrans(obj, data2Corr, transMat, refPlane, currentPlane)
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
            
            row = data2Corr.row;
            col = data2Corr.col;
            if ~isempty(idx2Corr)
                for j = 1:length(idx2Corr)
                                       
                    row = row + sign* transMat.rowTrans(idx2Corr(j));
                    col = col + sign* transMat.colTrans(idx2Corr(j));
                    
                end
                data2Corr.row = row;
                data2Corr.col = col;
                
            end
            
            corrData = data2Corr;
        
        end
        
        function [corrData] = applyRot(obj, data2Corr, corrMat, refPlane, currentPlane)
            
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
                    %Pad Z coordinate
                    data2C(:,1) = data2Corr.row;
                    data2C(:,2) = data2Corr.col;
                    data2C(:,3) = 0;
                    %remove center of mass (CM)
                    CM = mean(data2C);
                    data2C = data2C - CM;
                    
                    %Change the orientation of the data (should be [x;y;z]
                    %not [x y z]
                    data2C =  data2C';
                    
                    %Correction occur here
                    for j = 1:length(idx2Corr)
                        
                        if sign
                            rot = corrMat.rot{idx2Corr(j)}';
                        else
                            rot = corrMat.rot{idx2Corr(j)};
                        end
                           
                            data2C  = (rot*data2C);
                            data2Store = data2C';
                            corrData   = data2Store(:,1:2)+CM(1:2);
                    end
                else
                    corrData = data2Corr(:,1:2);
            
                end
  
        end
    
        function [zPos,inRange] = getZPosition(obj,val2Z,zCal,currentPlane,method)
            
            relZ = obj.calibrated.oRelZPos;

            zRange = zCal.fitZParam.zRange;
            zRange = zRange{currentPlane};
            zVec = zRange(1):1:zRange(2); %Here we assume accuracy >= 1nm

            switch method
                case 'poly'

                    fit = polyval(zCal.calib{currentPlane,1},zVec);

                case 'spline'
                    fit = ppval(zCal.calib{currentPlane,2},zVec);
            end

            %find the index of the value the closest to the particle
            %ellipticity
             [~,idx] = min(abs(fit-val2Z));

             zPos = zVec(idx)+ relZ(currentPlane)*1000;          
             inRange = and(val2Z>=zCal.fitZParam(1).ellipRange(1),...
                 val2Z<=zCal.fitZParam(1).ellipRange(2));
            
             if isempty(zPos)
                 disp('ouuups zpos is empty');
             end

        end
        
        function [zRange] = getZRange(obj,ellipRange,zCal,currentPlane,method)
            
            relZ = obj.calibrated.oRelZPos;
                       
            zVec = -2000:1:2000; %Here we assume accuracy >= 1nm
            
            switch method
                case 'poly'
                    
                    fit = polyval(zCal.calib{currentPlane,1},zVec);
                
                case 'spline'
                    fit = ppval(zCal.calib{currentPlane,2},zVec);
            end
            
            %find the index of the value the closest to the particle
            %ellipticity
            
             [~,idx1] = min(abs(fit-ellipRange(1)));
             [~,idx2] = min(abs(fit-ellipRange(2)));
             
             zPos1 = zVec(idx1)+ relZ(currentPlane)*1000;    
             zPos2 = zVec(idx2)+ relZ(currentPlane)*1000;
             
             zRange = [zPos1, zPos2];
             
        end
        
        function [data]  = resolveXYZ(obj,partData)
         
            pxSize = obj.info.pxSize;
            ellipRange = obj.ZCal.fitZParam.ellipRange;  

            idx2Keep = and(partData.ellip > ellipRange(1), partData.ellip < ellipRange(2));
            partData(~idx2Keep,:) = table(nan);

            row  = partData.row(3)*pxSize;
            col  = partData.col(3)*pxSize;
            z    = partData.z(3);
            data = table(row,col,z,'VariableNames',{'row','col','z'});
            %check how to perform averaging depending on the camera config
            [doAvg]  = obj.checkDoAverage(partData.ellip(3));

            if doAvg
                elliptRange = ellipRange(1):0.001:ellipRange(2);
                %we weigh the average later base on how much out of focus the
                %plane was.
                wRange1 = length(elliptRange(elliptRange<=1));
                wRange2 = length(elliptRange(elliptRange>=1));
                weight1 = linspace(1,5,wRange1);
                weight2 = linspace(5,1,wRange2);
                finalWeight = [weight1 weight2];
                ellipKept = partData.ellip(idx2Keep);
                idx = ellipKept;
                for k = 1 :length(ellipKept)

                    [~,idx(k)] = min(abs(elliptRange-ellipKept(k)));

                end

                weight = finalWeight(idx);
                %Weighed average
                row = sum(diag(partData.row(idx2Keep)* weight))/sum(weight) * pxSize;
                col = sum(diag(partData.col(idx2Keep)* weight))/sum(weight) * pxSize;
                z   = sum(diag(partData.z(idx2Keep)* weight))/sum(weight) * pxSize;
            end

            data.rowM = row;
            data.colM = col;
            data.zM = z;
            
        end
         
        function [data]  = resolveXYZInt(obj,partData,frameData)
           
            pxSize = obj.info.pxSize;
            ROIRad = ceil(obj.info.FWHM_px/2+1);
            
            planes  = partData(~isnan(partData.plane),:).plane;
            nPlanes = length(planes);
            bf = partData.plane(3);
            planePos = obj.calibrated.oRelZPos;
            if nPlanes==1
                z= 0;
            else
                %Get ROI XZ, YZ scaled to same pixel size
                [Mag] = Core.MPLocMovie.getZPhasorMag(partData,ROIRad,frameData);
                
                domain = 1:length([Mag.x]);
                data   = [Mag.x]+[Mag.y];
                guess.sig = 1.5*obj.info.FWHM_px;
                guess.mu  = bf;
                [Res,~] = SimpleFitting.gauss1D(data,domain,guess);
                
                z = Res(2);

            end
            %if the z position is out of bound we do not consider the data
            if or(z<min(domain),z>max(domain))
                z   = NaN;                           
                row = NaN;
                col = NaN;
                zM   = NaN;                           
                rowM = NaN;
                colM = NaN;
            else
                tmpZ = floor(z);
                fracZ = z-tmpZ;
                z = planePos(tmpZ)+fracZ*(planePos(tmpZ+1) - planePos(tmpZ));
                z = z*1000;
                                         
                row = partData.row(3)*pxSize;
                col = partData.col(3)*pxSize;
                zM = z;                      
                rowM = partData.row(3)*pxSize;
                colM = partData.col(3)*pxSize;
                
            end
            %store the data
            data = table(row,col,z,rowM,colM,zM,...
                   'VariableNames',{'row','col','z','rowM','colM','zM'});
            
        end
        
        function [method] = pickZFitMethod(obj)
            
            names = fieldnames(obj.ZCal.zAccuracy);
            nMethod = numel(names);
            
            if nMethod == 1
                method = names{1};
            else
                for i = 1: nMethod
                  currentMethod = names{i};
                  currentAccuracy = obj.ZCal.zAccuracy.(names{i}).BestFocus;
                  
                  %Here accuracy should be small (high accuracy mean small
                  %number)
                  if i==1
                    finalMethod = currentMethod;
                  elseif and(i>1, or(currentAccuracy > finalMethod,currentAccuracy==0))
                      
                      %finalMethod stay the same
                  
                  elseif and(i>1, and(currentAccuracy <  finalMethod,currentAccuracy>0))
                  
                      finalMethod = currentMethod;              
                       
                  end
                                                    
                end
                method = finalMethod;
            end
            
            
        end     
    end
    
    methods (Static)
        function [Mag] = getZPhasorMag(partData,ROIRad,volIm)

        %Possible improvement : Translate the coordinate of the best
        %focus into the otherplanes to extract the exact value where
        %the PSF should be taken    
        imSize = size(volIm);
        pos = [round(nanmean(partData.row)),round(nanmean(partData.col))];

        ROIs = Misc.getROIs(pos,ROIRad,imSize(1:2));

        ROI = volIm(ROIs(1):ROIs(2),ROIs(3):ROIs(4),:);

        Mag = struct('x',zeros(1,size(ROI,3)),'y',zeros(1,size(ROI,3)));
        for i =1:size(ROI,3)
            [~,~,~,Mag(i).x,Mag(i).y] = Localization.phasor(ROI(:,:,i));
        end

        end 

    end
end