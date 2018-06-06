classdef Movie <  handle
    %MOVIE will hold most of the information about the data and is the
    %object with which the user will mainly interact.
    
    %   Detailed explanation goes here
    
   properties (SetAccess = 'private')
       
       raw
       cal
       status = 'none';
       info
       calibrated
       superResCal
 
    end
    
    methods
        function obj = Movie(raw, info, cal, calibrated, superResCal)
            %MOVIE Construct an instance of this class
            %   We allow the user to create this object with various number
            %   of input allowing therefore to restart the analysis at any
            %   steps in the process.
            obj.raw = raw;
            %Give a status depending on input.
            switch nargin
                case 1 
                    obj.status = 'raw';
                case 2
                    obj.info = info;
                    obj.status = 'raw';
                case 3
                    obj.info = info;
                    obj.cal = cal;
                    obj.status = 'rawCal';
                case 4
                     obj.cal = cal;
                     obj.info = info;
                     obj.calibrated = calibrated;
                     obj.status = 'calibrated';
                case 5
                     obj.cal = cal;
                     obj.info = info;
                     obj.calibrated = calibrated;
                     obj.superResCal = superResCal;
                     obj.status = 'SRcalibrated';
                otherwise 
                    error('Unexpected number of argument');
            end
        end
        
        function  set.raw(obj,raw)
            [checkRes] = obj.checkRaw(raw);
            %When every check/modification is performed, raw is finally
            %store into the movie object
            if checkRes == false
                answer = questdlg('Path of raw file seems fine but something is wrong with the info, should we re-extract it?',...
                    'Question to user','Yes','No','Yes');
                switch answer
                    case 'Yes'
                        [~, movInfo, ~ ] = Load.Movie.ome.getInfo(raw.path);
                        raw.info = movInfo;
                    case 'No'
                        error('Something is wrong with raw file, request aborted.')
                end
                
            end
            obj.raw = raw;
        end
        
        function set.info(obj,info)
            %TO DO write the needed assertion.
            [checkRes] = obj.checkInfo(info);
            if checkRes
                obj.info = info;
            else
                error('Something is wrong with the format of the info provided by the user')
            end
        end
        
        function set.cal(obj,cal)
            [checkRes] = obj.checkCal(cal);
            if checkRes == false
                answer = questdlg('Path of calibration seems fine but something is wrong with the info, should we re-extract it?',...
                    'Question to user','Yes','No','Yes');
                switch answer
                    case 'Yes'             
                        [cal, movInfo] = mpSetup.cali.calculate(cal.path, false);
                        cal.info = movInfo;
                        cal.file = cal;
                    case 'No' 
                        error('Something is wrong with calibration, request aborted.')
                end
            end
            obj.cal = cal;
        end
        
        function set.status(obj,status)
            %TO DO Check with Rafa how 
            [checkRes] = obj.checkStatus(status);
            if checkRes
            obj.status = status;
            else
                error('There is a mismatch between the current experimental status and the content of the movie object, please check.')
            end 
        end

        function set.calibrated(obj,calibrated)
            %TO DO write the needed assertion.
            obj.calibrated = calibrated;
        end
        
        function set.superResCal(obj,superResCal)
            %TO DO write the needed assertion.
            obj.superResCal = superResCal;
        end
        
        function [checkRes] = checkRaw(~,raw)
            checkRes = false;
            %check that raw is not empty
            assert(~isempty(raw), 'Raw data should not be empty, a path should at least be provided');
            %Extract the fields names
            fields = fieldnames(raw);
            %Check that the first fields is path
            assert(and(isfield(raw,'path'),~isempty(raw.path)),...
                    'no filepath provided in Raw');
            %Check that the file path exist
            assert(exist(raw.path,'file') == 2,...
             'File specified in raw cannot be found');
            if (length(fields)==1)
            else %if 2 fields check that the format of raw.info is correct 
                if(or(isfield(raw,'info'),and(~isstruct(raw.info),...
                        or(isempty(raw.info),...
                        numel(fieldnames(raw.info)) ~= 6))))
                else
                    checkRes = true;
                end
            end
        end
    
        function [checkRes] = checkInfo(~,info)
            checkRes = true;
            if(isempty(info))
                warning('You did not provide info to the movie')
            else
                if ~isstruct(info)
                    checkRes = false;
                end
            end
        end
        
        function [checkRes] = checkCal(~,cal)
            checkRes = false;
            %check that raw is not empty
            assert(~isempty(cal), 'Calibration data should not be empty'); 
            %Extract the fields names
            fields = fieldnames(cal);
            %Check that the first fields is path
            assert(and(isfield(cal,'path'),~isempty(cal.cal)),...
                    'no filepath provided in Cal');
            %Check that the file path exist
            assert(exist(cal.path,'file') == 2,...
             'File specified in cal.path cannot be found')
         
            if (length(fields)<3)
            else
                if(or(and(isfield(cal,'info'),isfield(cal,'file')),...
                        or(and(~isstruct(cal.info),...
                        or(isempty(cal.info),numel(fieldnames(cal.info)) ~= 6)),...
                        and(~isstruct(cal.file),...
                        or(isempty(cal.file),numel(fieldnames(cal.file)) ~= 9)))))
                else
                    checkRes = true;
                end
            end
        end
         
        function [checkRes] = checkStatus(obj,status)
            checkRes = false;
            switch status
                case 'none'
                    if(isprop(obj,'raw'))
                    else
                        checkRes = true;
                    end
                case 'raw'
                    if(isprop(obj,'raw'))
                    [checkRes] = obj.checkRaw(obj.raw);
                    end
                case 'rawCal'
                    [checkRes] = obj.checkStatus('raw');
                    if(isprop(obj,'cal'))
                        [checkRes2] = obj.checkCal(obj.cal);
                        checkRes = checkRes*checkRes2;                        
                    end
                case 'Calibrated'
                    [checkRes] = obj.checkStatus('rawCal');
                        %add checkCalibrated
                   % checkRes = checkRes*checkRes2;                        
                   
                case 'SRCalibrated'
                    [checkRes] = obj.checkStatus('Calibrated');
                        %addCheckSRCalibrated
                        %checkRes = checkRes*checkRes2;                        
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    methods (Access = private)
        
        
        function updateStatus(obj)
            
            if(isempty(obj.raw))
                obj.status = 'none';
            elseif (isempty(obj.cal))
                obj.status = 'raw';
            elseif (isempty(obj.calibrated))
                obj.status = 'rawCalc';
            elseif (isempty(obj.superResCal))
                obj.status = 'calibrated';
            else
                obj.status = 'SRCalibrated';
            end

        end
        
    end
end

