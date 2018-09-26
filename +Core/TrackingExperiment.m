classdef TrackingExperiment < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        path
        trackMovies
        cal2D
        ZCal
        SRCal
        traces3D
        RMSD
        
    end
    
    methods
        function obj = TrackingExperiment(folder2Data,cal2D,SRCalPath,zCalPath)
            %TrackingExperiment Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = folder2Data;
            obj.cal2D = cal2D;
            obj.ZCal = zCalPath;
            obj.SRCal = SRCalPath;
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
                
                obj.SRCal.cal = cal;
                obj.SRCal.path = SRCal;
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
                
                obj.ZCal.cal = cal;
                obj.ZCal.path = zCal;
            end
        end
        
        function retrieveMovies(obj)
            %we get the zCalibration directory
            folder2Mov = dir(obj.path);
            folder2Mov = folder2Mov(cell2mat({folder2Mov.isdir}));
            %loop through the content of the directory
            for i = 3:size(folder2Mov,1)
                %Check if the directory
                currDir = dir([folder2Mov(i).folder filesep folder2Mov(i).name]);
                idx = contains({currDir.name}, 'ome.tif');
                if ~all(idx==0)
                    %we extract z motor position to check if the movie
                    %is indeed a zCalibration (expect zStack)
                    tmp = Core.MPTrackingMovie([folder2Mov(i).folder filesep folder2Mov(i).name], obj.cal2D, obj.SRCal.path, obj.ZCal.path);
                    obj.trackMovies.(['mov' num2str(i-2)]) = tmp;
                    
                else
                    
                    warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ome.Tif and is therefore ignored']);
                    
                end
                
            end
            disp('=======> DONE ! <========')
        end
        
        function retrieveTrackData(obj,detectParam, trackParam,val2Use)
            %Checking user input
            assert(nargin==4, 'retrieveZCalData expects 2 inputs, 1)detection Parameters, tracking parameter');
            assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
            assert(and(isfield(trackParam,'euDistXY'),isfield(trackParam,'euDistZ')),...
                'Tracking parameter is expected to be a struct with two field "euDistPXY" and "euDistZ"')
            
            %Extraction of Data
            nfields = numel(fieldnames(obj.trackMovies));
            allTraces = [];
            for i = 1: nfields
                
                disp(['Retrieving data from tracking file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                currentTrackMov = obj.trackMovies.(['mov' num2str(i)]);
                if i == 1
                    %Ask user for info about the setup for detection
                    currentTrackMov.giveInfo;
                    
                else
                    %get the info about the setup stored into the first
                    %object
                    currentTrackMov.info = obj.trackMovies.(['mov' num2str(1)]).getInfo;
                    
                end
                
                %Molecule detection
                currentTrackMov.findCandidatePos(detectParam);
                
                %SR fitting
                currentTrackMov.SRLocalizeCandidate;
                
                %apply SRCal

                currentTrackMov.applySRCal;
                
                %apply ZCal
                currentTrackMov.applyZCal;
                
                %Plane consolidation
                currentTrackMov.consolidatePlanes;
                
                %superResolve
                currentTrackMov.superResolve(val2Use);
                
                %tracking occurs here
                
                currentTrackMov.trackParticle(trackParam);
                
                currentTrackMov.getTraces;
                
                allTraces = [allTraces; traces];
                
            end
            
            obj.traces3D = allTraces;
            
            filename = [obj.path filesep 'traces3D.mat'];
            save(filename,allTraces);
            
            
            disp('=================> DONE <===================');
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
            obj.RMSD = MSD;
            obj.traces3D = traces;
                     
            filename = [obj.path filesep 'traces3D.mat'];
            save(filename,'traces');
            
            filename = [obj.path filesep 'RMSD-all.mat'];
            save(filename,'MSDmat');
            
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

