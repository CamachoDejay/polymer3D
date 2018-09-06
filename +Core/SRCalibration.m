classdef SRCalibration < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path
        SRCalMovies
        cal2D
        calib
        accuracy
    end
    
    methods
        
        function obj = SRCalibration(path2SRCal,cal2D)
            %zCalibration Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = path2SRCal;
            obj.cal2D = cal2D;
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
                        tmp = Core.SRCalMovie([folder2Mov(i).folder filesep folder2Mov(i).name], obj.cal2D);
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
            %Checking user input
            assert(nargin==3, 'retrieveSRCalData expects 3 inputs, 1)detection Parameters, fit z parameter, tracking parameter');
            assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
            assert(and(isstruct(trackParam),and(isfield(trackParam,'euDistPx'),isfield(trackParam,'ellip'))),'trackParam is expected to be a struct with 2 fields "euDistPx", the tolerated euclidian distance between consecutive frames, "ellip" ellipticity score to reach to be considered as frame consistent ([1 9])');
            
            %Extraction of SRData
            nfields = numel(fieldnames(obj.SRCalMovies));
            allData = [];
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
                
                %plane consolidation
                obj.SRCalMovies.(['SRCal' num2str(i)]).superResConsolidate;
                
%                 %frame consolidation/track In Z
%                 [traces,counter] = obj.SRCalMovies.(['SRCal' num2str(i)]).trackInZ(trackParam);
%                 
                %getting Calibration data
                [SRCalibData] = obj.SRCalMovies.(['SRCal' num2str(i)]).getSRCalData(trackParam);
                
                %Get data for every particle every planes together stored
                %in allData.
                
                    
                 allData = [allData ; SRCalibData ];

            end
            
             %Sort the data
            allData{1,3} = obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip{1,3};
            for j = 1: length(obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip)
                [~,ind] = sort(allData{j,1}(:,1));
                allData{j,1} = allData{j,1}(ind,:);
                
            end
            %Sort all data together
            [~,ind] = sort(allData{1,2}(:,1));
            allData{1,2} = allData{1,2}(ind,:);
            
            
            disp('=================> DONE <===================');
            %store the results of the data extraction
            obj.calib.trackParam = trackParam;
            obj.calib.data = allData;
        end
        
    end
end

