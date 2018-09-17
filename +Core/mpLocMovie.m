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
        
        function obj = MPLocMovie(raw, MPCal, SRCal, zCal)
            
            obj  = obj@Core.MPParticleMovie(raw,MPCal);
            
            switch nargin
                
                case 1
                    error('MPCal is required to create MPLocMovie')
                case 2
                    obj.SRCal = [];
                    obj.ZCal  = [];
                case 3
                    obj.SRCal = SRCal;
                    obj.ZCal = [];
                case 4
                    obj.SRCal = SRCal;
                    obj.ZCal = zCal;
                otherwise
                    error('Too many input arguments');
            end
        end
        
        function set.SRCal(obj,SRCal)
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
            
        end
        
        function set.ZCal(obj,zCal)
            
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
        end
        
        function applyCorr(obj,rot,refPlane)
            
            %apply SRCal
            obj.applySRCal(rot,refPlane);
            
            %transform ellipticity into Z
            obj.applyZCal;
            
        end
        
        function applySRCal(obj, rot, refPlane)
            assert(~isempty(obj.unCorrLocPos),'You need to find candidate and SR Localized them before applying corrections');
            assert(~isempty(obj.SRCal),'SR Calibration needed to correct the data');
            
            if nargin <2
                refPlane = 5;
            end
            if isempty(obj.corrLocPos)
                obj.corrLocPos = obj.localizedPos;
            end
            
            data = obj.unCorrLocPos;
            
            disp(['Applying SR calibration...']);
            for i = 1 : length(data)
                currData = data{i};
                currPlanes = unique(currData.plane);
                for j = 1 : length(currPlanes)
                    currentPlane = currPlanes(j);
                    data2Corr = currData(currData.plane==currentPlane,{'row','col','plane'});
                    
                    
                    if rot
                        corrMat = obj.SRCal.rot;
                        [corrData] = Core.SRCalMovie.applyRot(data2Corr, corrMat,refPlane);
                        
                    else
                        corrMat = obj.SRCal.trans;
                         [corrData] = Core.SRCalMovie.applyTrans(data2Corr,corrMat,refPlane);                    
                    end
                    
                    %we store the corrected data
                    obj.corrLocPos{i}(currData.plane==currentPlane,{'row','col','plane'}) = corrData;
                    
                end
               
            end
            obj.corrected.XY = true;
            disp('========> DONE ! <=========');

        end
        
        function applyZCal(obj)
            
            assert(~isempty(obj.unCorrLocPos),'You need to find candidate and SR Localized them before applying corrections');
            assert(~isempty(obj.ZCal),'Z Calibration needed to correct the data');
            
            if isempty(obj.corrLocPos)
                obj.corrLocPos = obj.localizedPos;
                warning('Z calibration is currently being applied on non-SRCorrected (X-Y) data');
            end
            
            data = obj.corrLocPos; 
            zCal = obj.ZCal;
            
            %we check which method is best:
            [method] = obj.pickZFitMethod;
            
            %Here we translate ellipticity into z position based on
            %calibration
            disp('Applying Z Calibration');
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
                
        function superResolve(obj, val2Use)
            disp('super resolving positions ... ');
            data2Resolve = obj.particles.List;
            
            SRList = cell(data2Resolve);
            for i = 1:length(data2Resolve)
                frameData = data2Resolve{i};
                for j = 1:length(frameData)
                    SRList{i}{j} = SRList{i}{j}(1,{'row','col','z'});
                    partData = frameData{j};
                    [row,col,z] = obj.resolveXYZ(partData(:,{'row','col','z','ellip'}),val2Use);
                    
                    SRList{i}{j}.row = row;
                    SRList{i}{j}.col = col;
                    SRList{i}{j}.z = z;
                    
                end
            end
                
            obj.particles.SRList = SRList;    
            disp('========> DONE ! <=========');
            
        end
                   
        function showCorrLoc(obj,frames)
             part = obj.particles.SRList;
            switch nargin
                case 1
                    frames = 1:length(part);
                case 2 
                    [frames] = obj.checkFrame(frames,obj.raw.maxFrame(1));
            end
           
            
            figure()
            hold on
            for i = 1:length(frames)
                cFrame = frames(i);
                currentFrame = part{cFrame};
                if ~isempty(length(part{cFrame}))
                    
                    for j = 1: length(currentFrame)
                        currentPart = currentFrame{j};
                        data = table2array(currentPart(1,{'row','col','z'}));
                        sizeMarker = 5;
                        scatter3(data(1),data(2),data(3),sizeMarker,data(3),'filled');
                       
                        
                    end
                end
                    
            end
            
            title('all Localization plotted');
            xlabel('x position in nm');
            ylabel('y position in nm');
            zlabel('z position in nm');
            
            
            hold off
            
        end
        
    end
    
    methods (Static)
        
                    
    end
    
    
    methods (Access = private)
        
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
    
        function [zPos,inRange] = getZPosition(obj,ellip,zCal,currentPlane,method)
            
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
             [~,idx] = min(abs(fit-ellip));
             
             zPos = zVec(idx)+ relZ(currentPlane)*1000;          
             inRange = and(ellip>=zCal.fitZParam(1).ellipRange(1),...
                 ellip<=zCal.fitZParam(1).ellipRange(2));
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
        
        function [row,col,z]  = resolveXYZ(obj,partData,val2Use)
            ellipRange = obj.ZCal.fitZParam.ellipRange;  
            pxSize = obj.info.pxSize;
            switch val2Use
                case 'bestFocus'
                   row = partData.row(3)*pxSize;
                   col = partData.col(3)*pxSize;
                   z = partData.z(3);
                case 'mean'
                    idx2Keep = and(partData.ellip > ellipRange(1), partData.ellip < ellipRange(2));
                    row = mean(partData.row(idx2Keep))*pxSize;
                    col = mean(partData.col(idx2Keep))*pxSize;
                    z   = mean(partData.z(idx2Keep));
                otherwise
                    error('Unknown value to use, only know "bestFocus" and "mean"');
            end
          
            
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
end