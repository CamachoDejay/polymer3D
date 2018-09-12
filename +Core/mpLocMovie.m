classdef MPLocMovie < Core.MPParticleMovie
    %MPLocMovie hold all the method and info necessary for localization,
    %expension of this object will be first centered around Tracking but
    %it will also be easy to extend to STORM/PALM.
    
    properties (SetAccess = 'protected')
        
        SRCal
        ZCal
        corrLocPos
    
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
                assert(istable(cal),...
                    'SR calibration is supposed to be a table');
                
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
            assert(~isempty(obj.localizedPos),'You need to find candidate and SR Localized them before applying corrections');
            assert(~isempty(obj.SRCal),'SR Calibration needed to correct the data');
            
            if nargin <2
                refPlane = 5;
            end
            if isempty(obj.corrLocPos)
                obj.corrLocPos = obj.localizedPos;
            end
            
            data = obj.localizedPos;
           corrMat = obj.SRCal;
            disp(['Applying SR calibration...']);
            for i = 1 : length(data)
                currData = data{i};
                currPlanes = unique(currData.plane);
                for j = 1 : length(currPlanes)
                    currentPlane = currPlanes(j);
                    data2Corr = currData(currData.plane==currentPlane,{'row','col'});
                    
                    [corrData] = obj.applyTrans(data2Corr,corrMat,refPlane,currentPlane);
                    
                    if rot 
                        [corrData] = obj.applyRot(corrData, corrMat,refPlane,currentPlane);
                    end
                    %we store the corrected data
                    obj.corrLocPos{i}(currData.plane==currentPlane,{'row','col'}) = corrData;
                    
                end
                %correctedData{i}(:,3) = data{i}(:,3);
            end
            disp('========> DONE ! <=========');
           % obj.corrLocPos = correctedData;
            
        end
        
        function applyZCal(obj)
            
            assert(~isempty(obj.localizedPos),'You need to find candidate and SR Localized them before applying corrections');
            assert(~isempty(obj.ZCal),'Z Calibration needed to correct the data');
            
            if isempty(obj.corrLocPos)
                obj.corrLocPos = obj.localizedPos;
                warning('Z calibration is currently being applied on non-SRCorrected (X-Y) data');
            end
            
            data = obj.localizedPos; 
            zCal = obj.ZCal;
            
            %Here we translate ellipticity into z position based on
            %calibration
            disp('Applying Z Calibration');
            for i = 1 : length(data)
                currData = data{i};
                nPos = size(currData,1);
                
                for j = 1 : nPos
                    
                    currentEllip = currData.ellip(j);
                    currentPlane = currData.plane(j);
                    [zPos] = getZPosition(obj,currentEllip,zCal,currentPlane);
                    obj.corrLocPos{i}.z(j) = zPos;
                    
                end
               
            end
            
            %Here we translate the ellipticity range into zRange for each
            %plane
            
            ellipRange = zCal.fitZParam.ellipRange;
            nPlanes = obj.calibrated.nPlanes;
            zRange = cell(nPlanes,1);
            for i = 1 : nPlanes
                zRange{i} = obj.getZRange(ellipRange,zCal,i);
            end
            disp('=======> DONE ! <========');
        end
        
        function [locPos] = getCorrLocPos(obj,frames)
             %Extract the position of the candidate of a given frame
            [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
            locPos = obj.corrLocPos{idx};
            
            if isempty(locPos)
                
                warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                
            end
        end
        
        function consolidatePlanes(obj,roiSize,frames)
            %Consolidation refers to connect molecules that were localized
            %at similar position in different plane on a single frame.
            assert(~isempty(obj.calibrated),'Data should be calibrated to consolidate');
            assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before consolidation');
            assert(~isempty(obj.localizedPos),'Localization needs to be performed before consolidation');
           
            %Check if some particles were saved already.
            [run, particle] = Core.MPParticleMovie.existParticles(obj.raw.movInfo.Path, '.mat');
            
            if run
                %Check the number of function input
                switch nargin
                    case 1
                        
                        frames = 1: obj.calibrated.nFrames;
                        roiSize = 6;
                        disp('Running consolidation on every frame with roi of 6 pixel');
                        
                    case 2
                        
                        frames = 1: obj.calibrated.nFrames;
                        disp('Running consolidation on every frame')
                        
                    case 3
                        
                        [frames] = Movie.checkFrame(frames,obj.raw.maxFrame(1));
                        assert(min(size(roiSize))==1,'RoiSize is expected to be a single number')
                        assert (isnumeric(roiSize),'RoiSize is expected to be a single number');
                        
                    otherwise
                        
                        error('Something wrong with input');
                        
                end
                
                nFrames = length(frames);
                %allocate for storage
                particleList = cell(1,obj.raw.maxFrame(1));
                nParticles = zeros(1,obj.raw.maxFrame(1));
                idx2TP = zeros(1,obj.raw.maxFrame(1));
                h = waitbar(0,'Consolidating candidate ...');
                
                %Consolidation occurs here
                for i = 1 : 1:nFrames
                    disp(['Consolidating frame ' num2str(i) ' / ' num2str(nFrames)]);
                    idx = frames(i);
                    %#1 Extract localized Position for specific frame
                    [fCandMet] = getCorrLocPos(obj,idx);
                    
                    if isempty(fCandMet)
                        
                        warning('Frame %d did not contain any localized positions',idx);
                        particleList{idx} = [];
                        %particleList{idx}{1} = nan(5);
                        nParticles(idx) = 0;
                        
                    else
                        %#2 Consolidate the position of the given frame
                        %across plane
                          %Calculate a focus metric (FM) combining ellipticity and GLRT FM.
                            [corrEllip, focusMetric] = Localization.calcFocusMetric(fCandMet.ellip,fCandMet.LRT);

                            %reformating to keep the same format as how the data is saved
                            %later
                            fCandMet.LRT = focusMetric;
                            fCandMet.Properties.VariableNames{5} = 'ellipLRT';
                            
%                             fCandMet = [fCandMet fCandMet(:,end)];
%                             fCandMet(:,6) = fCandMet(:,5);
%                             fCandMet(:,5) = focusMetric;
%                             fCandMet(:,4) = [];

                            focusMetric((1-corrEllip)>0.3) = NaN;
                            
                            %Plane Consolidation occur here
                            [part] = obj.planeConsolidation(fCandMet,focusMetric);

                            %we delete empty cells from the array
                            idx2Empty = cellfun(@isempty,part);
                            part(idx2Empty(:,1),:) = [];
                   
                            particleList{idx} = part;
                            nParticles(idx) = length(part);
                            
                            if ~isempty(part)
                                idx2TP(idx) = idx;
                            end
                        
                    end
                    waitbar(i/nFrames,h,['Consolidating candidate... ' num2str(i) '/' num2str(nFrames) ' done']);
                end
                close(h);
                
                %#3 Storing List
                particle.List       = particleList;
                particle.nParticles = nParticles;
                particle.tPoint     = nFrames;
                particle.idx2TP     = nonzeros(idx2TP);
                particle.roiSize    = roiSize;
                particle.Traces     = [];
                particle.nTraces    = [];
                
                fileName = sprintf('%s%sparticle.mat',obj.raw.movInfo.Path,'\');
                save(fileName,'particle');
            end
            %#4 Storing particles in the object
            obj.particles = particle;
        end
        
        function superResolve(obj)
            %TO BE PROGRAMMED
        end
                   
        function showCorrLoc(obj)
            part = obj.particles.List;
            
            figure()
            hold on
            for i = 1:length(part)
                currentFrame = part{i};
                if ~isempty(length(part{i}))
                    
                    for j = 1: length(currentFrame)
                        currentPart = currentFrame{j};
                        data = table2array(currentPart(3,{'row','col','z'}));
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
    
        function [zPos,inRange] = getZPosition(obj,ellip,zCal,currentPlane)
            
            relZ = obj.calibrated.oRelZPos;
                       
            zVec = -2000:1:2000; %Here we assume accuracy >= 1nm
                       
            fit = polyval(zCal.calib{currentPlane,1},zVec);
            %fit = ppval(zCal{particle(3,end),2},zVec);%spline
            
            %find the index of the value the closest to the particle
            %ellipticity
            
             [~,idx] = min(abs(fit-ellip));
             
             zPos = zVec(idx)+ relZ(currentPlane)*1000;          
             inRange = and(ellip>=zCal.fitZParam(1).ellipRange(1),...
                 ellip<= zCal.fitZParam(1).ellipRange(2));

        end
        
        function [zRange] = getZRange(obj,ellipRange,zCal,currentPlane)
            
            relZ = obj.calibrated.oRelZPos;
                       
            zVec = -2000:1:2000; %Here we assume accuracy >= 1nm
                       
            fit = polyval(zCal.calib{currentPlane,1},zVec);
            %fit = ppval(zCal{particle(3,end),2},zVec);%spline
            
            %find the index of the value the closest to the particle
            %ellipticity
            
             [~,idx1] = min(abs(fit-ellipRange(1)));
             [~,idx2] = min(abs(fit-ellipRange(2)));
             
             zPos1 = zVec(idx1)+ relZ(currentPlane)*1000;    
             zPos2 = zVec(idx2)+ relZ(currentPlane)*1000;
             
             zRange = [zPos1, zPos2];
             
        end
       
    end
end