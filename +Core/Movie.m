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
        function obj = Movie(raw, cal, info, calibrated, superResCal)
            %MOVIE Construct an instance of this class
            %   Detailed explanation goes here

            obj.raw = raw;
            obj.cal = cal;
      
            %Give a status depending on input.
            switch nargin
                case 1
                    obj.status = 'raw';
                case 2
                    obj.status = 'rawCal';
                case 3
                    obj.status = 'rawCal';
                    obj.info = info;
                case 4
                     obj.status = 'calibrated';
                     obj.info = info;
                     obj.calibrated = calibrated;
                case 5
                     obj.status = 'SRcalibrated';
                     obj.superResCal = superResCal;
                     obj.info = info;
                     obj.calibrated = calibrated;
                otherwise 
                    error('Unexpected number of argument');
            end
        end
        
        function  set.raw(obj,raw)
            %check that raw is not empty
            assert(~isempty(raw), 'Raw data should not be empty');
            %check that raw isa struct with no more than 2 fields
            assert(and(isstruct(raw),numel(fieldnames(raw))<3),...
                'raw should be a structure with no more than 2 fields');
            %Extract the fields names
            fields = fieldnames(raw);
            %Check that the first fields is path
            assert(and(isfield(raw,'path'),~isempty(raw.path)),...
                    'no filepath provided in Raw');
            %Check that the file path exist
            assert(exist(raw.path,'file') == 2,...
             'File specified in raw cannot be found')
            %if only one field extract and store movInfo into raw.info
            if (length(fields)==1)
                warning('Info missing in raw, extracting info')
                [~, movInfo, ~ ]= Load.Movie.ome.getInfo(raw.path);
                raw.info = movInfo;
            else %if 2 fields check that the format of raw.info is correct 
                %otherwise reextract it based on the raw.path
                assert(isfield(raw,'info'),'Unexpected field names in raw')
                if(and(~isstruct(raw.info),or(isempty(raw.info),...
                        numel(fieldnames(raw.info)) ~= 6)))
                  [~, movInfo, ~ ] = Load.Movie.ome.getInfo(raw.path);
                  warning('Info missing  or unexpected format in raw, re-extracting info based on the path')
                  raw.info = movInfo;
                end
            end
            %When every check/modification is performed, raw is finally
            %store into the movie object
            obj.raw = raw;
        end
        
        function set.cal(obj,cal)
            %check that raw is not empty
            assert(~isempty(cal), 'Calibration data should not be empty');
            %check that cal is a struct with no more than 2 fields
            assert(and(isstruct(cal),numel(fieldnames(cal))<4),...
                'Cal should be a structure with no more than 3 fields'); 
            %Extract the fields names
            fields = fieldnames(cal);
            %Check that the first fields is path
            assert(and(isfield(cal,'path'),~isempty(cal.(fields{1}))),...
                    'no filepath provided in Cal');
            %Check that the file path exist
            assert(exist(cal.path,'file') == 2,...
             'File specified in cal.path cannot be found')
         
            if (length(fields)<3)
                warning('Info and/or calibration missing in cal, extracting info and calibration file based on the path')
                [cal, movInfo] = mpSetup.cali.calculate(filePath, false);
                cal.info = movInfo;
                cal.file = cal;
            else
                assert(and(isfield(cal,'info'),isfield(cal,'file')),...
                'Something is wrong in the fieldnames in cal');
                if(or(and(~isstruct(cal.info),or(isempty(cal.info),...
                        numel(fieldnames(cal.info)) ~= 6)),and(~isstruct(cal.file),or(isempty(cal.file),numel(fieldnames(cal.file)) ~= 9))))
                    warning('something is wrong with cal.info and/or cal.file, re-extracting these based on the path')
                    [cal, movInfo] = mpSetup.cali.calculate(cal.path, false);
                    cal.info = movInfo;
                    cal.file = cal;
                end
            end
            obj.cal = cal;
        end
        
        function set.status(obj,status)
            %TO DO write check function to see if the data in the object
            %make sense with the current status.
            switch status
                case 'none'
                case 'raw'
                case 'rawCalc'
                case 'Calibrated'
            end
            obj.status = status;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

