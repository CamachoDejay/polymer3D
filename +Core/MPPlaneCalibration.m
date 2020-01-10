classdef MPPlaneCalibration < handle
    %PlaneCalibration is a class that will hold the information of several
    %MPCalMovie and be able to process these information to create a
    %calibration file from the data of all the movies
    
    properties (SetAccess = 'private')
        path
        ext
        MPCalibrations
        info
        allCal
        cal        
    end
    
    methods
        
        function obj = MPPlaneCalibration(path2MPCal,info)
            %zCalibration Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = path2MPCal.path;
            obj.ext  = path2MPCal.ext;
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
                folderPath = [folder2Mov(i).folder filesep folder2Mov(i).name];
                file2Analyze = Core.Movie.getFileInPath(folderPath,obj.ext);
               
                if ~isempty(file2Analyze)
                    
                    file.path = file2Analyze.folder;
                    file.ext  = obj.ext;
                    %we extract z motor position to check if the movie
                    %is indeed a zCalibration (expect zStack)
                    tmp = Core.MPCalibration(file, obj.info);
                    obj.MPCalibrations.(['MPCal' num2str(i-2)]) = tmp;

                else

                    warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ome.Tif and is therefore ignored']);

                end

            end
            disp('=======> DONE ! <========')
        end
        
        function calcIndivCal(obj)
            
            fieldN = fieldnames(obj.MPCalibrations);
            nfields = numel(fieldN);
      
            for i = 1: nfields
                
                if isfield(obj.info,'nChan')
                    nChan = obj.info.nChan;
                else
                    nChan = 4;
                end
                
                disp(['Retrieving data from MPCalibration file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                
                obj.MPCalibrations.(fieldN{i}).calc(nChan);
                currentCal = obj.MPCalibrations.(fieldN{i}).getCal;
                
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
            allFit = zeros([size(allData(1).file.fit) nFiles]);
            allNewOrder =  zeros([size(allData(1).file.neworder) nFiles]);
            allICorrF   =  zeros([size(allData(1).file.neworder) nFiles]);
            inFocus = allData(1).file.inFocus;
            inFocus = rmfield(inFocus,{'frame', 'zpos'});
            allRelZpos = zeros(nPlanes,nFiles);
            allZpos = zeros(nPlanes,nFiles);
            
            for i = 1:nFiles
                
                allROI(:,:,i) = allData(i).file.ROI;
                allFocusMet(:,:,i) = allData(i).file.focusMet;
                allFit(:,:,i) = allData(i).file.fit;
                allNewOrder(:,:,i) = allData(i).file.neworder; 
                allICorrF(:,:,i)  = allData(i).file.Icorrf;   
                tmp = cell2mat({allData(i).file.inFocus.zpos});
                allRelZpos(:,i) = tmp-tmp(1);
                allZpos(:,i) = tmp;
 
            end

            ROI = floor(mean(allROI,3));
            RelZPos = mean(allRelZpos,2);
            tmp = num2cell(RelZPos');
            [inFocus(1,:).relZPos] = tmp{:};
            tmp = num2cell(mean(allZpos,2)');
            [inFocus(1,:).zPos] = tmp{:};
            obj.cal.nFiles = nFiles;
            obj.cal.fullPath = obj.path;
            obj.cal.file = obj.allCal(1).file;
            obj.cal.file = rmfield(obj.cal.file,{'focusMet','fit','Zpos'});
            obj.cal.file.ROI = ROI;
            obj.cal.file.Icorrf = squeeze(mean(allICorrF,3));
            obj.cal.file.inFocus = inFocus;
            [~] = obj.determineCAMConfig;
            disp('================>DONE<====================');
            
            
            
        end
        
        function save(obj)
            cal2D = obj.getCal;
            filePath = obj.path;
            fileName = [filePath filesep '2DCal.mat'];
            
            save(fileName,'cal2D')
            disp('The calibration was succesfully saved');
        end
        
        function showCal(obj,idx)
            fields = fieldnames(obj.MPCalibrations);
            mov2Use = obj.MPCalibrations.(fields{idx});
            focusMet = mov2Use.cal.file.focusMet;
            fit = mov2Use.cal.file.fit(:,2:2:end);
            fitZ = mov2Use.cal.file.fit(:,1:2:end);
            ZPos = mov2Use.cal.file.Zpos;
            color = rand(8,3);
            
            FocusZ = {mov2Use.cal.file.inFocus.zpos};
            figure()
            hold on
            leg = cell(1,size(focusMet,2));
            height =  max(max(focusMet));
            y = 1:height;
            for i = 1 : size(focusMet,2)
                [~,idx] = max(fit(:,i));
                scatter(ZPos-mean(ZPos),focusMet(:,i),[],color(i,:),'filled')
                plot(fitZ(:,i)-mean(ZPos),fit(:,i),'Color', color(i,:),'LineWidth',2.5,'HandleVisibility','off')
               
                x = ones(1,length(y))*(FocusZ{i}-mean(ZPos));
                plot(x(:),y(:),'k--','HandleVisibility','off');
                
                leg{i} = ['Cam' num2str(obj.cal.file.inFocus(i).cam) ' - Plane' num2str(obj.cal.file.inFocus(i).ch)];
                
            end
            ylim([min(min(focusMet)), max(max(focusMet))]);
            xlim([round(min(ZPos-mean(ZPos))), round(max(ZPos-mean(ZPos)))]);
            title('Setup Plane Calibration');
            ylabel('Intensity (a.u.)')
            xlabel('z position (um)')
            legend(leg)
            hold off
    
        end
        
        function [camConfig] = determineCAMConfig(obj)
        relZPos = cell2mat({obj.cal.file.inFocus.relZPos});
        relZPos = relZPos(obj.cal.file.neworder);
        planeDist = abs(mean(diff(relZPos)))*1000;
        minDist = min(abs(relZPos(2:end)))*1000;
        if planeDist > 350
            camConfig = 'fullRange';
        elseif and(minDist>100, planeDist>200)
            camConfig = 'interleaved';
        elseif minDist<100
            camConfig = 'equal';
        else
            error('Something is wrong with your distance between planes')
        end
        
        obj.cal.camConfig = camConfig;

        end
        
        function offTarget(obj)
            newOrder = obj.cal.file.neworder;
            zPos = abs([obj.cal.file.inFocus.zPos]);
            zPos = zPos(newOrder);
            switch obj.cal.camConfig
                case 'fullRange'
                    distBetweenCamPlanes = abs(mean(diff(zPos(1:end/2))) + mean(diff(zPos(end/2+1:end))))/2;
                    target    = distBetweenCamPlanes;
                    distBetweenCam = abs(zPos(5)-zPos(4));
                    offTarget = distBetweenCam - target;
                   
                case 'interleaved'
                    distBetweenCamPlanes = abs(mean(diff(zPos(1:2:end))) + mean(diff(zPos(2:2:end))))/2;
                    target    = distBetweenCamPlanes/2;
                    distBetweenPlane = abs(diff(zPos));
                    offTarget1 = distBetweenPlane - target;
                    offTarget = mean(abs(offTarget1));

               
                case 'equal'
                    
                    distBetweenPlanes = abs(diff(zPos));
                    distBetweenPlanes = distBetweenPlanes(1:2:end);
                    offTarget = mean(distBetweenPlanes);
            end
           fprintf('The difference between the target and the current plane conformation \nis %d nm',round(offTarget*1000));

        end
        
        
        
     end
end