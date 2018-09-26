classdef imPlane < handle
    %IMPLANE hold single iamge plane information, e.g. a x-y frame of a z
    %stack.
    %   Detailed explanation goes here
    
    properties (Constant = true)
        posValidUnits = {'micron','nm'};
    end
    
    properties (SetAccess = 'private')
        data {mustBeNonnegative}
        numClass
        size {mustBeInteger}
        pixSizeNm {mustBeNonnegative}
        sizeNm
        position
        timeIdx
        timeStamp
    end
    
    methods
        function obj = imPlane(data)
            %IMPLANE Construct an instance of this class
            %   Detailed explanation goes here
            assert(ismatrix(data),'input data image must be a matrix of numeric non-negative values')
            obj.data = data;
            tmp = whos('data');
            obj.size = tmp.size;
            obj.numClass = tmp.class;
        end
        
        function setPixSizeNm(obj,value)
            obj.pixSizeNm = value;
            obj.sizeNm = obj.size.*value;
        end
        
        function setPosition(obj,x,y,z,units)
            
            % widths will be give by pix size in nm and image size!!!
            switch nargin
                case 2
                    error('Not enought iput arguments: x,y needed')
                case 3
%                     disp('x,y')
                    z = [];
                    units = 'none';
                case 4
%                     disp('x y z')
                    units = 'none';
                case 5
%                     disp('we can calculate widths')
                otherwise
                    error('unexpected')
                
            end
            mustBeNumeric(x)
            mustBeNonempty(x)
            mustBeNumeric(y)
            mustBeNonempty(y)
            mustBeNumeric(z)
                    
            obj.position.x = x;
            obj.position.y = y;
            obj.position.z = z;
            obj.position.xlims = nan(1,2);
            obj.position.ylims = nan(1,2);
            
            assert(any(strcmp(units,obj.posValidUnits)),'invalid units for postions, only "micron" and "nm" are accepted')
            
            switch units
                case 'micron'
                    % I need some conversion
                    obj.position.units = 'micron';
                    if ~isempty(obj.pixSizeNm)
                        % then I can calculate im lims
                        imWidthsMic = obj.sizeNm*1e-3;
                        dx = imWidthsMic(1)/2;
                        dy = imWidthsMic(2)/2;
                        obj.position.xlims = [x-dx, x+dx];
                        obj.position.ylims = [y-dy, y+dy];
                    end
                case 'nm'
                    % easier to handle
                    obj.position.units = 'nm';
                    if ~isempty(obj.pixSizeNm)
                        % then I can calculate im lims
                        imWidthsNm = obj.sizeNm;
                        dx = imWidthsNm(1)/2;
                        dy = imWidthsNm(2)/2;
                        obj.position.xlims = [x-dx, x+dx];
                        obj.position.ylims = [y-dy, y+dy];
                    end
                    
                otherwise
                    error('unexpected')
            end
            
            
            
        end
        
        
        function setTime(obj,idx,stamp)
            
            assert(and(idx>0, isinteger(idx)),'index must be a positive integer')
            
            switch nargin
                case 2
                    stamp = [];
                case 3
                    
                otherwise
                    error('Bad input')
            end
            
            obj.timeIdx = idx;
            obj.timeStamp = stamp;
            
        
        end
        
        
        % stters
        function set.pixSizeNm(obj,value)
            
            if value > 200
                warning('Are you sure the pixels size is correct, this is larger than ever expected')
            elseif value < 50
                warning('Are you sure the pixels size is correct, this is smaller than ever expected')
            end
            obj.pixSizeNm = value;
            
        end
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

