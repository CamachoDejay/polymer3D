classdef PlaneCalibration < handle
    %PlaneCalibration is a class that will hold the information of several
    %MPCalMovie and be able to process these information to create a
    %calibration file from the data of all the movies
    
    properties
        path
        MPCalibrations
        info
        allCal
        cal
        
    end
    
    methods
        
        function obj = PlaneCalibration(path2MPCal,info)
            %zCalibration Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = path2MPCal;            
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
        
        function retrieveMPCalibration(obj)
            %we get the MPCalibration directory
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
                        tmp = Core.MPCalibration([folder2Mov(i).folder filesep folder2Mov(i).name], obj.info);
                        tmp.calibrate;
                        [zStep, ~] = tmp.getZPosMotor;
                        %TODO: Check other motor position (we do not want
                        %any other movement here.
                        
                        if zStep ~= 0
                            %if it is we store
                            obj.MPCalibrations.(['MPCal' num2str(i-2)]) = tmp;
                        
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
            disp('=======> DONE ! <========')
        end
        
        function calcIndivCal(obj)
            
            nfields = numel(fieldnames(obj.MPCalibrations));
            allData = [];
            for i = 1: nfields
                
                if isfield(obj.info,'nChan')
                    nChan = obj.info.nChan;
                else
                    nChan = 4;
                end
                
                disp(['Retrieving data from MPCalibration file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                
                obj.MPCalibrations.(['MPCal' num2str(i)]).calc(nChan);
                currentCal = obj.MPCalibrations.(['MPCal' num2str(i)]).getCal;
                
                obj.allCal(i).file = currentCal.file;
                
            end
        end
        
        function calcCombinedCal(obj)
            allData = obj.allCal;
            
            
            
        end
        
        
    end
end