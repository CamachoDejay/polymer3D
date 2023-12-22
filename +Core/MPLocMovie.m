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
                if strcmp(obj.info.zMethod,'PSFE')
                    error('zMethod is selected is PSFE while no z calibration was provided')
                end
             
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
            
            if or(strcmp(zMethod,'Intensity'),strcmp(zMethod,'3DFit'))
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
           
            %Check if some particle were super resolved already:
            [run,SRList] = obj.existZResParticles(obj.info.runMethod,obj.raw.movInfo.Path,'.mat');
           
            if run
                data2Resolve = obj.particles.List;
                nPlanes = obj.calibrated.nPlanes;
                nParticles = sum(obj.particles.nParticles);
                pxSize = obj.info.pxSize;
                SRList = table(zeros(nParticles,1),...
                        zeros(nParticles,1), zeros(nParticles,1), zeros(nParticles,1),...
                        zeros(nParticles,1), zeros(nParticles,1), zeros(nParticles,1),...
                        zeros(nParticles,1), zeros(nParticles,1),zeros(nParticles,1),...
                        'VariableNames',...
                        {'row','col','z','rowM','colM','zM','adjR','intensity','SNR','t'});
                nFrames = length(data2Resolve);
                h = waitbar(0,'SuperResolving position...');

                for i = 1:nFrames
                  
                    frameData = data2Resolve{i};
                    frameData2Store = table(zeros(size(frameData)),...
                        zeros(size(frameData)),zeros(size(frameData)),zeros(size(frameData)),...
                        zeros(size(frameData)),zeros(size(frameData)),zeros(size(frameData)),...
                        zeros(size(frameData)),zeros(size(frameData)),zeros(size(frameData)),...
                        'VariableNames',...
                        {'row','col','z','rowM','colM','zM','adjR','intensity','SNR','t'});
                    
                    if strcmp(obj.info.zMethod,'PSFE')
                        for j = 1:length(frameData)
                            partData = frameData{j};
                            [data] = obj.resolveXYZ(partData(:,{'row','col','z','ellip','plane'}));
                            frameData2Store(j,{'row','col','z','rowM','colM','zM'}) = data;
                            frameData2Store.intensity(j) = partData.intensity(3);
                            frameData2Store.SNR(j) = partData.SNR(3);
                            frameData2Store.t(j) = i;
                        end
                        
                    else
                        
                        fData = obj.getFrame(i);
                        ROIRad = ceil(obj.info.FWHM_px/2+1);
                        
                        for j = 1:length(frameData)
                            partData = frameData{j};
                            
                            switch obj.info.zMethod
                                case 'Intensity'
                                    if nPlanes ==1
                                        row  = partData.row(3)*pxSize;
                                        col  = partData.col(3)*pxSize;
                                        z    = partData.z(3);
                                        rowM = partData.row(3)*pxSize;
                                        colM = partData.col(3)*pxSize;
                                        zM   = partData.z(3);
                                        adjR = 0; 
                                        data = table(row,col,z,rowM,colM,zM,adjR,...
                           'VariableNames',{'row','col','z','rowM','colM','zM','adjR'});

                                    else
                         
                                        partVolIm = obj.getPartVolIm(partData,ROIRad,fData);
                                        [data] = obj.resolveXYZInt(partData(:,{'row','col','z','ellip','plane'}),partVolIm);

                                    end

                                case '3DFit'
                                    if nPlanes ==1
                                        row  = partData.row(3)*pxSize;
                                        col  = partData.col(3)*pxSize;
                                        z    = partData.z(3);
                                        rowM = partData.row(3)*pxSize;
                                        colM = partData.col(3)*pxSize;
                                        zM   = partData.z(3);
                                        data = table(row,col,z,rowM,colM,zM,...
                           'VariableNames',{'row','col','z','rowM','colM','zM'});

                                    else
                                        [data] = obj.resolveXYZ3DFit(partData(:,{'row','col','z','ellip','plane'}),fData);
                                    end
                            end

                            frameData2Store(j,{'row','col','z','rowM','colM','zM','adjR'}) = data;
                            frameData2Store.intensity(j) = partData.intensity(3);
                            frameData2Store.SNR(j) = partData.SNR(3);
                            frameData2Store.t(j) = i;

                        end
                    end
                startIdx = find(SRList.row==0,1);   
                SRList(startIdx:startIdx+height(frameData2Store)-1,:) = frameData2Store;   
                waitbar(i/nFrames,h,['SuperResolving positions: frame ' num2str(i) '/' num2str(nFrames) ' done']);
                end
                profile('viewer')
                close(h);
                %clean up the list
                SRList(isnan(SRList.row),:) = [];
            end
            
            obj.particles.SRList = SRList;
            particle = obj.particles;
          
            %Save the data
            fileName = sprintf('%s%sparticle.mat',obj.raw.movInfo.Path,'\');
            profile('off')
            save(fileName,'particle');
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
         
        function [data] = resolveXYZInt(obj,partData,partVolIm)
          
            pxSize = obj.info.pxSize;
          
            bf = partData.plane(3);
            planePos = obj.calibrated.oRelZPos;
            
            %Get ROI XZ, YZ scaled to same pixel size
            [Mag] = Core.MPLocMovie.getZPhasorMag(partVolIm);

            domain = planePos;
            data   = [Mag.x]+[Mag.y];
            guess.sig = 2*obj.info.FWHM_px*pxSize/1000;
            guess.mu  = planePos(bf);
            guess.A   = max(data)-min(data);
           
%             [Res,fitData] = SimpleFitting.gauss1D(data,domain,guess);
%             RMS = sqrt(sum((data(:)-fitData(:)).^2)/length(data));
%             adjR = 1 - RMS.^2/var(data);
%           
            params = [guess.sig guess.mu guess.A min(data)];
         
            fun = @(x) SimpleFitting.minGauss1D(domain,data,x);
            opt = optimset('Display','off');
            % then we can do:
            [out, RMSD] = fminsearch(fun,params,opt);
            %normalize RMSD with mean data
            adjR = 1 - RMSD.^2/var(data);

            %[~, gaussFit] = SimpleFitting.minGauss1D(domain,data,out);
      
            z = out(2);
            %if the z position is out of bound we do not consider the data
            if or(z<min(domain),z>max(domain))
                z   = NaN;                           
                row = NaN;
                col = NaN;
                zM   = NaN;                           
                rowM = NaN;
                colM = NaN;
            else
                z = z*1000;
                row = partData.row(3)*pxSize;
                col = partData.col(3)*pxSize;
                zM = z;                      
                rowM = partData.row(3)*pxSize;
                colM = partData.col(3)*pxSize;

            end
           
            %store the data
            data = table(row,col,z,rowM,colM,zM,adjR,...
                   'VariableNames',{'row','col','z','rowM','colM','zM','adjR'});
            
   
        end
        
        function [data] = resolveXYZ3DFit(obj,partData,ROI)
             
            pxSize = obj.info.pxSize;
            bf = partData.plane(3);
            planePos = obj.calibrated.oRelZPos;

            x0 = mean([ROIs(3) ROIs(4)])*pxSize;
            y0 = mean([ROIs(1) ROIs(2)])*pxSize;
      
            z0 = planePos(bf)*1000;
            width.xy = 200;
            width.z  = 600;
            
            x = (ROIs(3):ROIs(4))*pxSize;
            y = (ROIs(1):ROIs(2))*pxSize;
            z = planePos*1000;
         
            [domX,domY,domZ] = meshgrid(x,y,z);
           
            dom(:,:,:,1) = domX;
            dom(:,:,:,2) = domY;
            dom(:,:,:,3) = domZ;

            %Multiple gaussian fitting occurs here
            [gPar,resnorm,res] = Localization.Gauss.MultipleGFit3D(ROI,x0,y0,z0,dom,1,width);
            z = gPar(:,7);

            if or(z<min(planePos*1000),z>max(planePos*1000))
                z   = NaN;                     
                row = NaN;
                col = NaN;
                zM   = NaN;                           
                rowM = NaN;
                colM = NaN;
            else
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
        function [run, SRList] = existZResParticles(runMethod,Path, ext)
            SRList = [];
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'particle')==true)
                        particle = load([file2Analyze(1).folder filesep 'particle.mat']);
                        particle = particle.particle;
                        if isfield(particle,'SRList')
                            if ~isempty(particle.SRList)
                                run = false;
                                SRList = particle.SRList;
                            else
                                run = true;
                            end
                        else
                            run = true;
                        end
                    else
                
                        run = true;
                        
                
                    end
                    
                case 'run'
                    
                     run = true;
                     
                     
            end
        end
        function [Mag] = getZPhasorMag(ROI)

            %Possible improvement : Translate the coordinate of the best
            %focus into the otherplanes to extract the exact value where
            %the PSF should be taken    
           
            Mag = struct('x',zeros(1,size(ROI,3)),'y',zeros(1,size(ROI,3)));
            for i =1:size(ROI,3)
                [~,~,~,Mag(i).x,Mag(i).y] = Localization.phasor(ROI(:,:,i));
            end

        end
        
        function [partVolIm] = getPartVolIm(partData,ROIRad,volIm)
            %Extract data from particle in the 8 planes
            imSize = size(volIm);
            
            pos = [round(nanmean(partData.row)),round(nanmean(partData.col))];

            ROIs = Misc.getROIs(pos,ROIRad,imSize(1:2));

            partVolIm = volIm(ROIs(1):ROIs(2),ROIs(3):ROIs(4),:);
            
        end

    end
end