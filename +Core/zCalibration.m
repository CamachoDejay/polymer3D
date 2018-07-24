classdef zCalibration < handle
    %zCalibration is a class that will hold the information of several
    %zCalMovie and be able to process these information to create a
    %calibration curve together with displaying etc...
    
    properties
        path
        zCalMovies
        cal2D
        calib
        traces3D
        zPosMotor
        zAccuracy
    end
    
    methods
        
        function obj = zCalibration(path2zCal,cal2D)
            %zCalibration Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = path2zCal;
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
        
        function retrieveZCalMov(obj)
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
                        tmp = Core.zCalMovie([folder2Mov(i).folder filesep folder2Mov(i).name], obj.cal2D);
                        [zStep, ~] = tmp.getZPosMotor;
                        %TODO: Check other motor position (we do not want
                        %any other movement here.
                        
                        if zStep ~= 0
                            %if it is we store
                            obj.zCalMovies.(['zCal' num2str(i-2)]) = tmp;
                        
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
        
        function retrieveZCalData(obj,detectParam, fitZParam, trackParam)
            %Checking user input
            assert(nargin==4, 'retrieveZCalData expects 3 inputs, 1)detection Parameters, fit z parameter, tracking parameter');
            assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
            assert(and(isstruct(fitZParam),and(isfield(fitZParam,'deg'),isfield(fitZParam,'ellipRange'))),'fitZParam is expected to be a struct with 2 fields "deg" for the degree of polynomial parameter, "ellipRange" holding the min and max value of accepted ellipticity. ');
            assert(and(isstruct(trackParam),and(isfield(trackParam,'euDistPx'),isfield(trackParam,'ellip'))),'trackParam is expected to be a struct with 2 fields "euDistPx", the tolerated euclidian distance between consecutive frames, "ellip" ellipticity score to reach to be considered as frame consistent ([1 9])');
            
            %Extraction of zData
            nfields = numel(fieldnames(obj.zCalMovies));
            allData = cell(8,3);
            for i = 1: nfields
                disp(['Retrieving data from zCal file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                if i == 1
                    %Ask user for info about the setup for detection
                    obj.zCalMovies.(['zCal' num2str(i)]).giveInfo;
                    
                else
                    %get the info about the setup stored into the first
                    %object
                    obj.zCalMovies.(['zCal' num2str(i)]).info = obj.zCalMovies.(['zCal' num2str(1)]).getInfo;
                    
                end
                %Molecule detection
                obj.zCalMovies.(['zCal' num2str(i)]).findCandidatePos(detectParam);
                
                %plane consolidation
                obj.zCalMovies.(['zCal' num2str(i)]).superResConsolidate;
                
                %frame consolidation/track In Z
                [traces,counter] = obj.zCalMovies.(['zCal' num2str(i)]).trackInZ(trackParam);
                
                %getting Calibration data
                [zCalData] = obj.zCalMovies.(['zCal' num2str(i)]).getCalData(traces,counter);
                
                %Synchronizing calibration data
                [~] = obj.zCalMovies.(['zCal' num2str(i)]).syncZCalData(zCalData,fitZParam);
                
                %Get data for every particle every planes together stored
                %in allData.
                for j = 1: length(obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip)
                    
                    allData{j,1} = [allData{j,1} ;obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip{j,1} ];
                    
                end
                
                allData{1,2} = [allData{1,2} ;obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip{1,2} ];
                
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
            obj.calib.fitZParam = fitZParam;
            obj.calib.trackParam = trackParam;
            obj.calib.data = allData;
        end
        
        function zCalibrate(obj)
            assert(~isempty(obj.calib),'No calibration data was found, please run retrieveZCalData before calibrating');
            assert(isfield(obj.calib,'data'),'No calibration data was found, please run retrieveZCalData before calibrating');
            assert(iscell(obj.calib.data),'Something is wrong with the calibration data stored, please run retrieveZCalData before calibrating');
            
            %Get the data from the object
            allData = obj.calib.data;
            %Fitting of the data occurs here
            zCalibration = obj.calZCalibration(allData);
            %storage to the object
            obj.calib.file = zCalibration;

        end
        
        function showZCalibration(obj)
            %Checking object before going further
            assert(~isempty(obj.calib),'No calibration data was found found, please run  "zCalibrate" before display');
            assert(isfield(obj.calib,'data'),'No calibration data was found, please run retrieveZCalData before calibrating');
            assert(iscell(obj.calib.data),'Something is wrong with the calibration data stored, please run retrieveZCalData before calibrating');
            assert(isfield(obj.calib,'file'),'No zCalibration found, please run zCalibrate before showZCalibration');
            assert(iscell(obj.calib.file),'No zCalibration found, please run zCalibrate before showZCalibration');
            
            %Plot #1
            relZ = obj.zCalMovies.('zCal1').calibrated.oRelZPos*1000;%in nm
            figure()
            for i = 1 : length(obj.calib.data)
                
                z =  obj.calib.data{i}(:,1)+relZ(i);
                ellip = obj.calib.data{i}(:,2);
                
                ellip1 = ellip(ellip>=1);
                z1 = z(ellip>=1);
                
                ellip2 = 1./ellip(ellip<=1);
                z2 = z(ellip<=1);
                yAx = 1:0.1:2;
                subplot(2,1,1)
                hold on
                scatter(z1,ellip1)
                plot(ones(length(yAx),1)*relZ(i),yAx)
                title('Elliptiticy Elongated in Y')
                xlabel('zPos (nm)')
                ylabel('Ellipticity (sigY/sigX)')
                ylim([1 2])
                hold off
                
                subplot(2,1,2)
                hold on
                scatter(z2,ellip2)
                plot(ones(length(yAx))*relZ(i),yAx)
                title('Elliptiticy Elongated in X')
                xlabel('zPos (nm)')
                ylabel('Ellipticity (sigX/sigY)')
                ylim([1 2])
                hold off
                
            end
            
            %Plot #2
            figure()
            scatter(obj.calib.data{1,2}(:,1),obj.calib.data{1,2}(:,2));
            xlabel('ZPosition')
            ylabel('Ellipticity')
            title('Ellipticity-Z curve for all the planes superimposed')
            
            %Plot #3
            figure()            
            for i = 1 : length(obj.calib.data)
                
                dataCurrentPlane = obj.calib.data{i};
                
                [binnedData] = Plotting.qBinning(dataCurrentPlane,length(dataCurrentPlane)/5);
                zVec = min(dataCurrentPlane(:,1)):max(dataCurrentPlane(:,1));
                %retrieving fit to display
                p = obj.calib.file{i};
                fit = polyval(p,zVec);
                %shifting according to the plane
                zVec = zVec + relZ(i) ;
                dataCurrentPlane(:,1) = dataCurrentPlane(:,1)+ relZ(i);
                binnedData(:,1) = binnedData(:,1) +relZ(i);
                
                subplot(1,2,1)
                hold on
                scatter(binnedData(:,1),binnedData(:,2))
                plot(zVec,fit)
                title('Binned data fitted with polynomial')
                
                subplot(1,2,2)
                hold on
                scatter(dataCurrentPlane(:,1),dataCurrentPlane(:,2))
                plot(zVec,fit)
                title('Full data fitted with polynomial')
                
            end
            
            
        end
        
        function evalAccuracy(obj)
            
            nfields = numel(fieldnames(obj.zCalMovies));
            zMotor = cell(nfields,1);
            trace3D  = cell(nfields,1);
            
            for i = 1: nfields
                
                disp(['Retrieving 3D traces ' num2str(i) ' / ' num2str(nfields) ' ...']);
                [~,zMotor{i}] = obj.zCalMovies.(['zCal' num2str(i)]).getZPosMotor;
                trace3D{i} = obj.zCalMovies.(['zCal' num2str(i)]).get3DTraces (obj.calib.file,obj.calib.fitZParam);
                
            end
            obj.traces3D = trace3D;
            obj.zPosMotor = zMotor;
            
            disp('Calculating accuracy and preparing plots')
            obj.showAccuracy(trace3D,zMotor);
            
            disp('=================> DONE <===================');
        end
        
        
    end
    
    methods (Access = private)
        
        function [zCalibration] = calZCalibration(~,zSyncCalData)
            %function that take the synchronized z-ellip Data and fit, for
            %each planes with a polynomial. It stores the coeff of the
            %polynomials
            zCalibration = cell(length(zSyncCalData),1);
            deg = zSyncCalData{1,3}(3);
            
            for i = 1: length(zSyncCalData)
                dataCurrentPlane = zSyncCalData{i};
                [binnedData] = Plotting.qBinning(dataCurrentPlane,length(dataCurrentPlane)/5);
                
                zVec = min(dataCurrentPlane(:,1)):max(dataCurrentPlane(:,1));
                
                p = polyfit(dataCurrentPlane(:,1),dataCurrentPlane(:,2),deg);
                fit = polyval(p,zVec);
                
                zCalibration{i} = p;
                
            end
            
        end
        
        function showAccuracy(obj,traces3D, zPosMotor)
       
            %Display the x-y-z traces
            assert(~isempty(obj.traces3D),'You need to get the traces before displaying them, use evalAccuracy to get the traces');
            trace = traces3D;
            motor =  zPosMotor;
            accuracyFocus = [];
            accuracyMean = [];
            %plot XYZ for every particles       
            nfields = numel(fieldnames(obj.zCalMovies));
            figure()
            for i = 1: nfields
                
                currentTrace = trace{i};
                currentMotor = motor{i};
                currentMotor = (currentMotor - currentMotor(1))*1000;
                npart = size(currentTrace,3);
                
                for j = 1:npart
                    data = currentTrace(:,:,j);
                    if ~all(data==0)

                        data = data(data(:,1)~=0,:);
                        accuracyF2plot = (data(:,3) - data(find(data(:,3),1,'first'),3)) - currentMotor(1:size(data,1));
                        accuracyM2plot  = (data(:,6) - data(find(data(:,3),1,'first'),6)) - currentMotor(1:size(data,1));
                        
                        accuracyMean = [accuracyMean std(accuracyM2plot)];
                        accuracyFocus = [accuracyFocus std(accuracyF2plot)];
                          
                        subplot(2,2,1)
                        hold on
                        scatter(1:size(data(:,3),1),data(:,3) - data(1,3));
                        plot(currentMotor,'-b');

                        title('Z position for different particle in best focus')
                        xlabel('Frame')
                        ylabel('Position(nm)')
                        hold off

                        subplot(2,2,2)
                        hold on
                        scatter(1:size(data,1),data(:,6) - data(1,6));
                        plot(currentMotor,'-b');
                        xlabel('Frame')
                        ylabel('Position(nm)')
                        title('Z position for different particle mean')
                        hold off

                        subplot(2,2,3)
                        hold on
                        plot(nonzeros(accuracyF2plot))
                        xlabel('Frame')
                        ylabel('tracked Distance - motor movement')
                        title(' Z position for different particle in best focus')
                        hold off

                        subplot(2,2,4)
                        hold on
                        plot(nonzeros(accuracyM2plot))
                        xlabel('Frame')
                        ylabel('tracked Distance - motor movement')
                        title('Z position for different particle mean')
                        hold off
                    end
                end
            end
            
            disp(['Accuracy in Z for best focus is ' num2str(mean(accuracyFocus))]);
            disp(['Accuracy in Z for mean  ' num2str(mean(accuracyMean))]);
            
            obj.zAccuracy.BF = mean(accuracyFocus);
            obj.zAccuracy.M  = mean(accuracyMean);
        end
        
    end
end

