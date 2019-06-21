classdef TrackingExperiment < handle
    %This object purpose is to serve as a small app to analyze many
    %tracking data acquired in same conditions, combine results and make
    %some statistics out of the data.
    
    properties
        
        path
        ext
        trackMovies
        cal2D
        info
        ZCal
        SRCal
        traces3D
        MSD
        
    end
    
    methods
        function obj = TrackingExperiment(folder2Data,cal2D,info,SRCalPath,zCalPath)
            %TrackingExperiment Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = folder2Data.path;
            obj.ext  = folder2Data.ext;
            obj.cal2D = cal2D;
            obj.info = info;
            obj.ZCal = zCalPath;
            obj.SRCal = SRCalPath;
            
        end
        
        %Set function
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
            if isempty(cal2D)
                obj.cal2D = [];
            else
                assert(ischar(cal2D), 'Path should be given as a string');
                assert(isfolder(cal2D), 'The path given is not a folder, ZCalibration expect a folder. In the folder it is expected to find separate folder for each zCalMovie.')

                [file2Analyze] = Core.Movie.getFileInPath(cal2D,'2DCal.mat');

                if isempty(file2Analyze)
                    error('No 2D calibration file found in the given folder');
                else
                    fileName = [file2Analyze.folder filesep file2Analyze.name];
                    cal = load(fileName);
                    field = fieldnames(cal);
                    cal = cal.(field{1});
                    assert(and(isstruct(cal), and(isfield(cal,'camConfig'),isfield(cal,'file'))),...
                        '2D calibration is supposed to be a struct with 4 fields');
                    obj.cal2D = cal;

                end
            end

        end
        
        function set.SRCal(obj,SRCal)
            if isempty(SRCal)
                obj.SRCal.cal = SRCal;
                obj.SRCal.path = SRCal;
            else
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

                    obj.SRCal.cal = cal;
                    obj.SRCal.path = SRCal;
                end
            end
            
        end
        
        function set.ZCal(obj,zCal)
            if isempty(zCal)
                obj.ZCal.cal = zCal;
                obj.ZCal.path = zCal;
            else
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

                    obj.ZCal.cal = cal;
                    obj.ZCal.path = zCal;
                end
            end    
            end
        
        %get 3D traces
        
        function [traces3D] = getTraces3D(obj)
            traces3D = obj.traces3D;
        end
        
        %Extract movie from path
        
        function retrieveMovies(obj)
            %we get the zCalibration directory
            folder2Mov = dir(obj.path);
            folder2Mov = folder2Mov(cell2mat({folder2Mov.isdir}));
            %loop through the content of the directory
            for i = 3:size(folder2Mov,1)
                %Check if the directory
                folderPath = [folder2Mov(i).folder filesep folder2Mov(i).name];
                file2Analyze = Core.Movie.getFileInPath(folderPath,obj.ext);
               
                if ~isempty(file2Analyze)
                    
                    file.path = file2Analyze.folder;
                    file.ext  = obj.ext;
                    tmp = Core.MPTrackingMovie(file , obj.cal2D,obj.info, obj.SRCal.path, obj.ZCal.path);
                    tmp.calibrate;
                    obj.trackMovies.(['mov' num2str(i-2)]) = tmp;
                    
                else
                    
                    warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ome.Tif and is therefore ignored']);
                    
                end
                
            end
            disp('=======> DONE ! <========')
        end
        
        %Calculate tracking traces
        
        function retrieveTrackData(obj,detectParam, trackParam,val2Use)
            %Checking user input
            assert(nargin==4, 'retrieveZCalData expects 2 inputs, 1)detection Parameters, tracking parameter');
            assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
            assert(and(isfield(trackParam,'radius'),isfield(trackParam,'memory')),...
                'Tracking parameter is expected to be a struct with two field "radius" and "memory"')
            fieldsN = fieldnames(obj.trackMovies);
            %Extraction of Data
            nfields = numel(fieldsN);
            allTraces = [];
            for i = 1: nfields
                
                disp(['Retrieving data from tracking file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                currentTrackMov = obj.trackMovies.(fieldsN{i});
                if i == 1
                    %Ask user for info about the setup for detection
                    currentTrackMov.giveInfo;
                    
                else
                    %get the info about the setup stored into the first
                    %object
                    currentTrackMov.info = obj.trackMovies.(fieldsN{1}).getInfo;
                    
                end
                
                %Molecule detection
                currentTrackMov.findCandidatePos(detectParam);
                
                %SR fitting
                currentTrackMov.SRLocalizeCandidate;
                refPlane = round(currentTrackMov.calibrated.nPlanes/2);
                rot = true;
                %apply SRCal
                currentTrackMov.applySRCal(rot,refPlane);
                
                %apply ZCal
                currentTrackMov.applyZCal;
                
                %Plane consolidation
                currentTrackMov.consolidatePlanes;
                
                %superResolve
                currentTrackMov.superResolve;
                
                %tracking occurs here
                currentTrackMov.trackParticle(trackParam);
                
                [traces] = currentTrackMov.getTraces;
                fileN = cell(length(traces),1);
                fileN(:,1) = {i};
                
                [xStep,xMotor] = currentTrackMov.getXPosMotor;
                [yStep,yMotor] = currentTrackMov.getYPosMotor;
                [zSt,zMotor]   = currentTrackMov.getZPosMotor;
                
                colMot = cell(length(traces),1);
                colMot(:,1) = {xMotor};
                colStep = cell(length(traces),1);
                colStep(:,1) = {xStep};
                
                rowMot = cell(length(traces),1);
                rowMot(:,1) = {yMotor};
                rowStep = cell(length(traces),1);
                rowStep(:,1) = {yStep};
                
                zMot = cell(length(traces),1);
                zMot(:,1) = {zMotor};
                zStep = cell(length(traces),1);
                zStep(:,1) = {zSt};
                
                allTraces = [allTraces; traces(:), fileN,colStep,colMot,rowStep,rowMot,zStep,zMot ];
                
            end
            
            obj.traces3D = allTraces;
            
%             filename = [obj.path filesep 'traces3D.mat'];
%             save(filename,'allTraces');
            
            
            disp('=================> DONE <===================');
        end
        
        %Plotting for individual movies
         
        function showLoc(obj,idx)
             fieldsN = fieldnames(obj.trackMovies);
             maxIdx = length(fieldsN);
             assert(idx <= maxIdx,['Requested index to Movie is too large, only ' num2str(maxIdx) ' movies']);
             
             currentMov = obj.trackMovies.(fieldsN{idx});
             
             currentMov.showCorrLoc;
        end
        
        function showTraces(obj,idx)
             fieldsN = fieldnames(obj.trackMovies);
             maxIdx = length(fieldsN);
             assert(idx <= maxIdx,['Requested index to Movie is too large, only ' num2str(maxIdx) ' movies']);
             
             currentMov = obj.trackMovies.(fieldsN{idx});
             
             currentMov.showTraces;
        end
        
        function evalAccuracy(obj,dim,idx)
            
            fieldsN = fieldnames(obj.trackMovies);
            maxIdx = length(fieldsN);
            assert(idx <= maxIdx,['Requested index to Movie is too large, only ' num2str(maxIdx) ' movies']);
            
            currentMov = obj.trackMovies.(fieldsN{idx});
            
            currentMov.evalAccuracy(dim);
            
        end
        
        function [int,SNR] = getAvgIntensity(obj)
            assert(~isempty(obj.traces3D),'You need to extract 3D traces before extracting average intensity');
            traces = obj.traces3D;
            nTraces = length(traces);
            int = zeros(nTraces,1);
            SNR = zeros(nTraces,1);
            for i = 1: length(traces)
                currentTrace = traces{i};
                int(i) = mean(currentTrace.intensity);
                SNR(i) = mean(currentTrace.SNR);
                
            end
            
            int = mean(int);
            SNR = mean(SNR);
           
            
        end
        
        function [MSD,traces] = getRMSD(obj,dimension)
            
            assert(~isempty(obj.traces3D),'You need to extract 3D traces before extracting RMSD');
            
            traces = obj.traces3D;
            
            switch nargin
                case 1
                    
                    dimension = '3D';
                    
                case 2
                    
                otherwise
                    
                    error('too many input arguments');
                    
            end
            MSD = cell(traces);
            MSDmat = zeros(obj.trackMovies.('mov1').raw.movInfo.maxFrame(1),length(traces));
            for i = 1 : length(traces)
                
                currentTrace = traces{i};
                coord = [currentTrace.col,currentTrace.row,currentTrace.z];
                [MSD,~] = Core.MPTrackingMovie.calcMeanSqrD(coord,dimension);
                currentTrace.RMSD = zeros(size(currentTrace.row,1),1);
                currentTrace.RMSD(1:end-1) = MSD;
                traces{i} = currentTrace;
                MSD{i} = MSD;
                MSDmat(1:length(MSD),i) = MSD(:);
                
            end 
            obj.MSD = MSD;
            obj.traces3D = traces;
                     
            
            
            
        end
        
        function saveData(obj)
            
            trackRes = struct; 
            disp('Saving Data');
            
            if ~isempty(obj.traces3D)
                
                trackData = obj.traces3D;
                MSDs = obj.MSD;
                    
                if ~isempty(MSDs)
                    trackRes.MSD = MSDs;
                end
                
                trackRes.traces = trackData; 
                trackRes.info = obj.info;
                trackRes.path = obj.path;
                filename = [obj.path filesep 'trackResults.mat'];
                save(filename,'trackRes');
                disp('Data were succesfully saved');
    
            else
                
                warning('No Data was saved because no traces or MSD could be found, please make sure you ran the analysis first');
            
            end   
        end
        
        function showRMSD(obj)
            MSD = obj.RMSD;
            
            figure()
            hold on
            
            for i = 1:length(MSD)
                currentMSD = MSD{i};
                plot(currentMSD)
            end
                       
        end
    end
end

