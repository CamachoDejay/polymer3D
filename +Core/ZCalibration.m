classdef ZCalibration < handle
    %zCalibration is a class that will hold the information of several
    %zCalMovie and be able to process these information to create a
    %calibration curve together with displaying etc...
    
    properties
        path
        zCalMovies
        cal2D
        info
        calib
        traces3D
        zPosMotor
        zAccuracy
        
    end
    
    methods
        
        function obj = ZCalibration(path2zCal,cal2D,info)
            %zCalibration Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = path2zCal;
            obj.cal2D = cal2D;
            obj.info = info;
            %We prepare zAccuracy
            obj.zAccuracy.spline = table(0,0,'VariableNames',{'BestFocus','Mean'});
            obj.zAccuracy.poly   = table(0,0,'VariableNames',{'BestFocus','Mean'});
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
            
            obj.cal2D = cal2D;
            
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
                        tmp = Core.MPZCalMovie([folder2Mov(i).folder filesep folder2Mov(i).name], obj.cal2D,obj.info);
                        tmp.calibrate;
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
            disp('=======> DONE ! <========')
        end
        
        function retrieveZCalData(obj,detectParam, fitZParam, trackParam)
            %Checking user input
            nPlanes = obj.zCalMovies.(['zCal' num2str(1)]).calibrated.nPlanes;
            assert(nargin==4, 'retrieveZCalData expects 3 inputs, 1)detection Parameters, fit z parameter, tracking parameter');
            assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
            assert(and(isstruct(fitZParam),and(isfield(fitZParam,'deg'),isfield(fitZParam,'ellipRangeCal'))),'fitZParam is expected to be a struct with 2 fields "deg" for the degree of polynomial parameter, "ellipRange" holding the min and max value of accepted ellipticity. ');
            assert(and(isstruct(trackParam),and(isfield(trackParam,'euDistPx'),isfield(trackParam,'commonPlanes'))),'trackParam is expected to be a struct with 2 fields "euDistPx", the tolerated euclidian distance between consecutive frames, "commonPlanes" nplane consistent to be considered as frame consistent ([1 4])');
            obj.calib.fitZParam.ellipRangeCal = fitZParam.ellipRangeCal;
            %Extraction of zData
            nfields = numel(fieldnames(obj.zCalMovies));
            allData = cell(nPlanes,3);
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
                
                %SR fitting
                obj.zCalMovies.(['zCal' num2str(i)]).SRLocalizeCandidate;
                
                %plane consolidation
                obj.zCalMovies.(['zCal' num2str(i)]).consolidatePlanes;
                
                %frame consolidation/track In Z
                [traces,counter] = obj.zCalMovies.(['zCal' num2str(i)]).trackInZ(trackParam);
                
                %getting Calibration data
                [zCalData] = obj.zCalMovies.(['zCal' num2str(i)]).getCalData(traces,counter);
                
                %Synchronizing calibration data
                [~] = obj.zCalMovies.(['zCal' num2str(i)]).syncZCalData(zCalData,fitZParam);
                
                %Get data for every particle every planes together stored
                %in allData.
                for j = 1: size(obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip,1)
                    
                    allData{j,1} = [allData{j,1} ;
                    obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip{j,1} ];
                    
                end
                
                allData{1,2} = [allData{1,2} ;obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip{1,2} ];
                
            end
            
             %Sort the data
            allData{1,3} = obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip{1,3};
            for j = 1: size(obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip,1)
                if ~isempty(allData{j,1})
                [~,ind] = sort(allData{j,1}.z);
                allData{j,1} = allData{j,1}(ind,:);
                end
            end
            %Sort all data together
            [~,ind] = sort(allData{1,2}.z);
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
        
        function showZCalibration(obj,method)
            %Checking object before going further
            assert(~isempty(obj.calib),'No calibration data was found found, please run  "zCalibrate" before display');
            assert(isfield(obj.calib,'data'),'No calibration data was found, please run retrieveZCalData before calibrating');
            assert(iscell(obj.calib.data),'Something is wrong with the calibration data stored, please run retrieveZCalData before calibrating');
            assert(isfield(obj.calib,'file'),'No zCalibration found, please run zCalibrate before showZCalibration');
            assert(iscell(obj.calib.file),'No zCalibration found, please run zCalibrate before showZCalibration');
            
            %Plot #1
            relZ = obj.zCalMovies.('zCal1').calibrated.oRelZPos*1000;%in nm
            zRange = obj.calib.fitZParam.zRange;
            ellipRange = obj.calib.fitZParam.ellipRangeCal;
            figure()
            for i = 1 : size(obj.calib.data,1)
                
                z =  obj.calib.data{i}.z+relZ(i);
                ellip = obj.calib.data{i}.ellip;
                
                ellip1 = ellip(ellip>=1);
                z1 = z(ellip>=1);
                
                ellip2 = 1./ellip(ellip<=1);
                z2 = z(ellip<=1);
                yAx = 1:0.1:2;
                subplot(2,1,1)
                hold on
                scatter(z1,ellip1)
                plot(ones(length(yAx),1)*relZ(i),yAx,'--k')
                title('Elliptiticy Elongated in Y')
                xlabel('zPos (nm)')
                ylabel('Ellipticity (sigY/sigX)')
                ylim([1 2])
                xlim([-4000,1000]);
                hold off
                
                subplot(2,1,2)
                hold on
                scatter(z2,ellip2)
                plot(ones(length(yAx))*relZ(i),yAx,'--k')
                title('Elliptiticy Elongated in X')
                xlabel('zPos (nm)')
                ylabel('Ellipticity (sigX/sigY)')
                ylim([1 2])
                xlim([-4000,1000]);
                hold off
                
            end
            
            %Plot #2
            figure()
            hold on
          
            for i = 1 : size(obj.calib.data,1)
                dataCurrPlane = obj.calib.data{i};
                scatter(dataCurrPlane.z, dataCurrPlane.ellip,15,'filled',...
                    'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4,'DisplayName',['Plane ' num2str(i) ' - ' num2str(relZ(i))])
               
            %scatter(obj.calib.data{1,2}(:,1),obj.calib.data{1,2}(:,2));
            end
            xlabel('ZPosition')
            ylabel('Ellipticity')
            legend(gca,'show')
            title('Ellipticity-Z curve for all the planes superimposed')
            hold off
            
            %Plot #3
            figure()   
            
            for i = 1 : size(obj.calib.data,1)
                
                dataCurrPlane = obj.calib.data{i};
                idx2Keep = and(dataCurrPlane.ellip>ellipRange(1),...
                    dataCurrPlane.ellip<ellipRange(2));
                
                [Res] = Core.ZCalibration.findConsecVal(idx2Keep);
 
                dataCurrPlane = dataCurrPlane(Res,:);
                [binnedData] = Plotting.qBinning([dataCurrPlane.z,...
                    dataCurrPlane.ellip],length(dataCurrPlane.z)/7);
                
                zVec = zRange{i}(1):zRange{i}(2);
               
                switch method
                    case 'poly'
                         p = obj.calib.file{i};
                         fit = polyval(p,zVec);
                    case 'spline'
                         p   = obj.calib.file{i,2};
                         fit = ppval(p,zVec);
                
                    otherwise
                        error('unknown method requested by the user');
                end
                %retrieving fit to display
               
                %shifting according to the plane
                zVec = zVec + relZ(i) ;
                dataCurrPlane.z = dataCurrPlane.z+ relZ(i);
                binnedData(:,1) = binnedData(:,1) +relZ(i);
                zVec = zVec(and(fit<ellipRange(2),fit>ellipRange(1)));
                fit = fit(and(fit<ellipRange(2),fit>ellipRange(1)));
                subplot(1,2,1)
                hold on
                markerSize = 10;
                scatter(binnedData(:,1),binnedData(:,2))
                plot(zVec,fit,'r')
                title('Binned data fitted with Spline')
                
                subplot(1,2,2)
                hold on
                scatter(dataCurrPlane.z,dataCurrPlane.ellip,markerSize,'filled')
                plot(zVec,fit,'r','LineWidth',2)
                title('Full data fitted with polynomial')
                
            end
        end
        
        function evalAccuracy(obj,fittingType)
            obj.calib.fitZParam.fittingType = fittingType;
            nfields = numel(fieldnames(obj.zCalMovies));
            %Allocate memory
            zMotor = cell(nfields,1);
            trace3D  = cell(nfields,1);
            %Loop through the movies
            for i = 1: nfields
                %info to user
                disp(['Retrieving 3D traces ' num2str(i) ' / ' num2str(nfields) ' ...']);
                [~,zMotor{i}] = obj.zCalMovies.(['zCal' num2str(i)]).getZPosMotor;
                trace3D{i} = obj.zCalMovies.(['zCal' num2str(i)]).get3DTraces (obj.calib.file,obj.calib.fitZParam,fittingType);
                
            end
                
            obj.traces3D = trace3D;
            obj.zPosMotor = zMotor;
            
            disp('Calculating accuracy and preparing plots')
            obj.showAccuracy(trace3D,zMotor);
            
            disp('=================> DONE <===================');
        end
        
        function [cal] = getCal(obj)
            
            cal = obj.calib;
            
        end
        
        function save(obj)
            cal.fitZParam = obj.calib.fitZParam;
            cal.fitZParam(1).type = 'polynomial';
            cal.fitZParam(2).type = 'spline';
            cal.calib = obj.calib.file;
            cal.zAccuracy = obj.zAccuracy;
            
            fileName = sprintf('%s%szCalibration.mat',obj.path,'\');
            save(fileName,'cal');
            disp('The calibration was succesfully saved');

        end
    end
    methods (Static)
         function [Res] = findConsecVal(bool)
            %determine the longuest section of consecutive true values
            %In this case, bool is mostly generated by using the conditions
            %on ellipticity range but it could be used in other cases.
            %!!Assume the longuest section is more or less centered!!
                    %Divide the data in 2 part
                    i = 1;
                    longestSec = [];
                    while i<length(bool)
                        if bool(i) == 1
                            startIdx = i;
                            endIdx = i;
                        
                        run = true;
                        while run
                            i = i+1;
                            if i == length(bool)
                                break;
                            end
                            if bool(i) == 1
                                endIdx = i;
                            else
                                run = false;
                            end
                        end
                        tmpSec = startIdx:endIdx;
                        
                        if length(tmpSec)> length(longestSec)
                            longestSec = tmpSec;
                            
                        end
                        
                        if endIdx==length(bool)
                            i = endIdx;
                        else
                            i = endIdx+1;
                        end
                        else
                            i = i+1;
                        end
                    end
                    Res = longestSec;                    
                    test = bool(Res);
                    if all(test==1)
                    else
                        disp('Something went wrong with findConsecVal');
                    end

        end
    end
    methods (Access = private)
        
        function [zCalibration] = calZCalibration(obj,zSyncCalData)
            %function that take the synchronized z-ellip Data and fit, for
            %each planes with a polynomial. It stores the coeff of the
            %polynomials
            ellipRange = obj.calib.fitZParam.ellipRangeCal;
            zCalibration = cell(length(zSyncCalData),2);
            deg = zSyncCalData{1,3}(3);
            disp('Starting fitting ...');
            for i = 1: size(zSyncCalData,1)
                disp(['Fitting of plane ' num2str(i)]);
                dataCurrPlane = zSyncCalData{i};
                idx2Keep = and(dataCurrPlane.ellip>ellipRange(1),...
                    dataCurrPlane.ellip<ellipRange(2));
                
                [Res] = Core.ZCalibration.findConsecVal(idx2Keep);
 
                dataCurrPlane = dataCurrPlane(Res,:);
                [binnedData] = Plotting.qBinning([dataCurrPlane.z,...
                    dataCurrPlane.ellip],length(dataCurrPlane.z)/7);
                
                %zVec = min(dataCurrentPlane.z):max(dataCurrentPlane.z);
                
                p = polyfit(dataCurrPlane.z,dataCurrPlane.ellip,deg);
                p2 = spline(binnedData(:,1),binnedData(:,2));
                %fit = polyval(p,zVec);
                
                zCalibration{i} = p;
                zCalibration{i,2} = p2;
                zRange = [min(dataCurrPlane.z), max(dataCurrPlane.z)];
                obj.calib.fitZParam.zRange{i} = zRange;
            end
            disp('=================> DONE <===================');
        end
        
        function showAccuracy(obj,traces3D, zPosMotor)
       
         %Display the x-y-z traces
            assert(~isempty(obj.traces3D),'You need to get the traces before displaying them, use evalAccuracy to get the traces');
            trace = traces3D;
            motor =  zPosMotor;
            accuracyZFocus = [];
            accuracyZMean = [];
            accuracyXFocus = [];
            accuracyYFocus = [];
            %plot XYZ for every particles       
            nfields = numel(fieldnames(obj.zCalMovies));
            zPlot =[];
            figure()
            for i = 1: nfields
                zStep = diff(motor{i});
                currentTrace = trace{i};
                currentMotor = motor{i};
                currentMotor = (currentMotor - currentMotor(1))*1000;
                npart = size(currentTrace,3);
                
                for j = 1:npart
                    data = currentTrace(:,:,j);
                    if ~all(data==0)
                        frameVec = find(data(:,1)~=0);
                        %frameVec(data(:,1)~=0,:) = [];
                        data = data(data(:,1)~=0,:);
                        bFit = mean(data(:,3)-zStep(1)*1000.*frameVec(:));
                        data2SubZ = zStep(1)*1000*frameVec+bFit;
                        
                        accuracyF2plot = (data(:,3)- data2SubZ(:));
                        accuracyX2plot = data(:,1) - mean(data(:,1));
                        accuracyY2plot = data(:,2) - mean(data(:,2));
                        
                        zPlot = [zPlot; data(:,3), accuracyF2plot];
                        
                        bFit = mean(data(:,6)-zStep(1)*1000.*frameVec(:));
                        data2SubZavg = zStep(1)*1000*frameVec+bFit;
                        accuracyM2plot  = (data(:,6)-data2SubZavg(:));
                        
                        accuracyZMean = [accuracyZMean mean(abs(accuracyM2plot))];
                        accuracyZFocus = [accuracyZFocus mean(abs(accuracyF2plot))];
                        
                       
                        accuracyXFocus = [accuracyXFocus mean(abs(accuracyX2plot))];
                        accuracyYFocus = [accuracyXFocus mean(abs(accuracyY2plot))];
                        
                        subplot(2,2,1)
                        hold on
                        scatter(1:size(data(:,3),1),data(:,3) );
                        plot(data2SubZ,'-b');

                        title('Z position for different particle in best focus')
                        xlabel('Frame')
                        ylabel('Position(nm)')
                        hold off

                        subplot(2,2,2)
                        hold on
                        scatter(1:size(data,1),data(:,6));
                        plot(data2SubZavg,'-b');
                        xlabel('Frame')
                        ylabel('Position(nm)')
                        title('Z position for different particle mean')
                        hold off

                        subplot(2,2,3)
                        hold on
                        plot(accuracyF2plot)
                        xlabel('Frame')
                        ylabel('tracked Distance - motor movement')
                        title(' Z position for different particle in best focus')
                        hold off

                        subplot(2,2,4)
                        hold on
                        plot(accuracyM2plot)
                        xlabel('Frame')
                        ylabel('tracked Distance - motor movement')
                        title('Z position for different particle mean')
                        hold off
                    end
                end
            end
            
            disp(['Accuracy in Z for best focus is ' num2str(nanmean(accuracyZFocus))]);
            disp(['Accuracy in X for best focus is ' num2str(nanmean(accuracyXFocus))]);
            disp(['Accuracy in Y for best focus is ' num2str(nanmean(accuracyYFocus))]);
            disp(['Accuracy in Z for mean  ' num2str(nanmean(accuracyZMean))]);
            
            switch obj.calib.fitZParam.fittingType
                case 'poly'
                    obj.zAccuracy.poly.('BestFocus') = mean(accuracyZFocus);
                    obj.zAccuracy.poly.Mean      = mean(accuracyZMean);
                case 'spline'
                    obj.zAccuracy.spline.BestFocus = mean(accuracyZFocus);
                    obj.zAccuracy.spline.Mean      = mean(accuracyZMean);
                otherwise
                    error('Unknown fitting type used, only currently known are "poly" and "spline"');
            end
            
            [~,idx] = sort(zPlot(:,1));
            zPlot = zPlot(idx,:);
            
           
            
            minBin = min(zPlot(:,1)):100:max(zPlot(:,1));
            maxBin = minBin + minBin(2)-minBin(1);
            midBin = (maxBin+minBin)/2;

            meanEllip = zeros(1,length(minBin));
            stdEllip   = zeros(1,length(minBin));

            for i=1:length(minBin)

            meanEllip(i) = 0;
            stdEllip(i) = mean(abs(zPlot(and(zPlot(:,1)>=minBin(i),zPlot(:,1)>=minBin(i)),2)));
            end


            figure
            Plotting.shadedErrorBar(midBin, meanEllip, stdEllip,'lineprops','-r','transparent',1)
            
            xlabel('Z position (nm)')
            ylabel('Mean error')
            
            
            
        end
        
       
    end
end

