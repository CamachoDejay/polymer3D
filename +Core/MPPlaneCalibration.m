classdef MPPlaneCalibration < handle
    %PlaneCalibration is a class that will hold the information of several
    %MPCalMovie and be able to process these information to create a
    %calibration file from the data of all the movies
    
    properties (SetAccess = 'private')
        path
        MPCalibrations
        info
        allCal
        cal        
    end
    
    methods
        
        function obj = MPPlaneCalibration(path2MPCal,info)
            %zCalibration Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = path2MPCal;            
            obj.info = info;
      
        end
        
        function [cal] = getCal(obj)
            cal = obj.cal;
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
        
        function retrieveMovies(obj)
            disp('Retrieving movies from indicated folder...')
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
                        obj.MPCalibrations.(['MPCal' num2str(i-2)]) = tmp;
 
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
            disp('Combining data from different calibration...');
            allData = obj.allCal;
            nFiles = length(allData);
            nPlanes = length(allData(1).file.neworder);
            allROI = zeros([size(allData(1).file.ROI) nFiles]);
            allFocusMet =  zeros([size(allData(1).file.focusMet) nFiles]);
            allFit = allFocusMet;
            allNewOrder =  zeros([size(allData(1).file.neworder) nFiles]);
            allICorrF   =  zeros([size(allData(1).file.neworder) nFiles]);
            inFocus = allData(1).file.inFocus;
            inFocus = rmfield(inFocus,{'frame', 'zpos'});
            allRelZpos = zeros(nPlanes-1,nFiles);
            
            for i = 1:nFiles
                
                allROI(:,:,i) = allData(i).file.ROI;
                allFocusMet(:,:,i) = allData(i).file.focusMet;
                allFit(:,:,i) = allData(i).file.fit;
                allNewOrder(:,:,i) = allData(i).file.neworder; 
                allICorrF(:,:,i)  = allData(i).file.Icorrf;
                
                tmp = cell2mat({allData(i).file.inFocus.zpos});
                allRelZPos(:,i) = tmp-tmp(1);
 
            end

            ROI = round(mean(allROI,3));
            RelZPos = mean(allRelZPos,2);
            test = num2cell(RelZPos');
            [inFocus(1,:).relZPos] = test{:};
            obj.cal.nFiles = nFiles;
            obj.cal.fullPath = obj.path;
            obj.cal.file = obj.allCal(1).file;
            obj.cal.file = rmfield(obj.cal.file,{'focusMet','fit','Zpos'});
            obj.cal.file.ROI = ROI;
            obj.cal.file.inFocus = inFocus;
            disp('================>DONE<====================');
        end
     end
end