classdef DDMMovie < Core.MPMovie
properties (SetAccess = 'private')
    
    IsEnoughRam
    IsCudaDevice

end
properties (SetAccess ='public')
    DDMOutput
    AllFrames
    DDMInfo
end
methods 
    
%%  Instantiation
        function obj = DDMMovie(raw, MPCal, info,DDMInfo)
             obj  = obj@Core.MPMovie(raw,MPCal,info); 

             assert(~isempty(DDMInfo), 'No MovieInfo provided, check parameter input'); 
             
             obj.DDMInfo = DDMInfo;

             try
                 gpuDevice;
                 obj.IsCudaDevice = true;         
             catch
                 obj.IsCudaDevice = false;
             end

        end
        


%%  Analysis




      

        
        

        
        function  [RadialValueInQSpace,ValidIndeces]=Get3DGrid(obj,sizes,critangle)
            
        % Prepares the grid for 3d averaging, with proper values for the
        % wavectors.
            try
              zpixel = abs(mean(diff(obj.calibrated.oRelZPos)));     
            catch
              zpixel = 1;   
            end

            for i=1:length(sizes)-1
                x{i} = 2*pi*(-round(sizes(i)/2):1:round(sizes(i)/2)-1)*1/obj.DDMInfo.PixelSize*1/sizes(i);
            end
            if sizes(end)~=1
            z = 2*pi*(-round(sizes(end)/2):1:round(sizes(end)/2)-1)*(1/zpixel)*1/sizes(end);
            else
            z = 0;   
            end
                
                
            [X,Y,Z] = meshgrid(x{2},x{1},z);
            RadialValueInQSpace = sqrt(X.^2 + Y.^2 + Z.^2);
            ValidIndeces =  Z./sqrt(X.^2 + Y.^2 + Z.^2)<abs(cosd(critangle));
        end
        
        
        

        
        
        
        function  LoadAllFrames(obj,driftCorr)
        % Loads all data. Stops if it almost runs out of memory ,- when the
        % amount of memory left is less then 2 GB, in order to leave some
        % of it for other operations. 
        
            obj.AllFrames = cell(1,obj.DDMInfo.nFrames);
            for i=1:obj.DDMInfo.nFrames 
                if rem(i,50)==0
                    disp([num2str(i) ' frames loaded']);
                end
                user = memory;
                if  user.MemAvailableAllArrays>1e+09
                    try
                    ZeroFrame= single(obj.getFrame(i));
                    catch
                    ZeroFrame= single(obj.getFrame(i).Cam1);
                    end
                    close all

                    
                    for k=1:length(size(ZeroFrame))
                       if mod(size(ZeroFrame,1),2) == 1 
                           ZeroFrame = ZeroFrame(2:end,:,:);
                       end
                       ZeroFrame = shiftdim(ZeroFrame,1);        
                    end
                    
                    obj.AllFrames{i}= ZeroFrame;  
                    
                else
                    obj.DDMInfo.nFrames = i-1;
                    warning(['Not enough RAM to load all of the '  num2str(obj.DDM.nFrames) 'frames. Loaded only ' num2str(i-1)] )
                    break;
                end
            end 
            %TODO: Fix data
            correlationInfo.corrSz = 100; %in px. Radius of the ROI used for correlation
            %correlation function
            correlationInfo.driftPeriod = 1; %in Frame, Number of frame that are averaged
            %for driftCalculation ==> 1 mean that drift is calculated for each frame
            scalingFactor = 1;%Used for interpolation in sub-pixel Drift correction 
            %objects of interest
            
            
            if driftCorr
                [corrData,~] = PreProcess.CorrelationDrift(obj.AllFrames,scalingFactor,correlationInfo);

                obj.AllFrames = corrData;
            end
            
        end
        
        
       function AvgFFT =  CalculateDelta(obj,AvgFFT,FrameSize,dt,ROI)
            for t=1:obj.DDMInfo.nFrames-dt-1
                if obj.IsCudaDevice==1
                    FrameDelta =gpuArray(obj.AllFrames{t+dt}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2))-obj.AllFrames{t}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2))); 
                else
                    FrameDelta =(obj.AllFrames{t+dt}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2))-obj.AllFrames{t}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2)));  
                end

                FrameDelta = abs( fftshift(fftn(FrameDelta-mean(FrameDelta,'all'),[FrameSize(1),FrameSize(2),FrameSize(3)]))).^2;                  
                AvgFFT = AvgFFT + FrameDelta;                                            
            end
        end 
        
       
        
            
        
        function  DDMOutput = main(obj,varargin)
         
            %Parse inputs conditionally
            p = inputParser;  
            addOptional(p, 'ROI',  [1, size(obj.AllFrames{1},1), size(obj.AllFrames{1},1);  %Default ROI as image size
                                    1, size(obj.AllFrames{1},2), size(obj.AllFrames{1},2);
                                    1 ,size(obj.AllFrames{1},3) ,size(obj.AllFrames{1},3)]);
            addOptional(p, 'Padsize',zeros(1, 3));
            addOptional(p, 'NumBins', 200);
            addOptional(p, 'CriticalAngle',0);
            parse(p,varargin{:});
            FrameSize = p.Results.ROI(:,3)'+p.Results.Padsize;

            DDMOutput = [];
            
            
            %Generate struct grid for n-D averaging

            for dt=1:obj.DDMInfo.FramesToAnalyze
                
                if obj.IsCudaDevice==1
                    AvgFFT = zeros(FrameSize(1),FrameSize(2),FrameSize(3),'gpuArray');
                else
                    AvgFFT = zeros(FrameSize(1),FrameSize(2),FrameSize(3));
                end
                    
                AvgFFT = obj.CalculateDelta(AvgFFT,FrameSize,dt, p.Results.ROI)./(obj.DDMInfo.nFrames-dt); 
                
                RadiallyAveragedDDMSignal=  obj.AverageRadialy3D(AvgFFT,p.Results.NumBins,FrameSize, p.Results.CriticalAngle);
                DDMOutput(:,1) =[NaN ; RadiallyAveragedDDMSignal(:,1)];
                DDMOutput(:,end+1) = [dt ; RadiallyAveragedDDMSignal(:,2)];
                
                disp([ num2str(dt) '/' num2str(obj.DDMInfo.FramesToAnalyze) ' Frames done'] );
                
            end
            close all
            DDMOutput =  obj.ConvertOutput(DDMOutput);
            obj.DDMOutput = DDMOutput;
        end
        
        
        
        function AverageDDMValueAtR = AverageRadialy3D(obj, AvgFFT,NumBins,FrameSize,critangle)
            
        % Averages the scattering function in 3d radially. Returns the
        % average values in function of time and q-vector.

                
                [RadialValueInQSpace, ValidRange]=obj.Get3DGrid(FrameSize,critangle);
                BinSize = max(RadialValueInQSpace,[],'all')/NumBins;
                FoundRadii = 0;
                R = 0;
                next = 1;
                AverageDDMValueAtR = zeros(NumBins,2);
                AvgFFT = gather(AvgFFT);
                while ~isempty(FoundRadii)
                    R = R+BinSize;
                    
                    FoundRadii = find(R-BinSize<=RadialValueInQSpace & RadialValueInQSpace<=R &  ValidRange );
                    
                    AverageDDMValueAtR(next,:) = [R ,  mean(AvgFFT(FoundRadii),'all')];
                    next = next+1;
                end
            

        end

     
end 
methods(Static)
    

        function DDMOutput = ConvertOutput(DDMOutList)
            DDMOutput = table(0,{[]},{[]},'VariableNames',{'QVector','Time','DDMSignalValue'});
            warning('off','all')
            for i=2:size(DDMOutList,1)
                DDMOutput.QVector(i-1) = DDMOutList(i,1);
                DDMOutput.Time(i-1) = {DDMOutList(1,2:end)};
                DDMOutput.DDMSignalValue(i-1) ={DDMOutList(i,2:end)};
    
            end
            warning('on','all')
        end
     



end  
end