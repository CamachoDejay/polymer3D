classdef FitDDMMovie < handle
%% Class to fit DDM output data
    
    properties (SetAccess ='public')
    FitResults
    DDMDataFit
    Model
    end
    
methods   
    
        function obj = FitDDMMovie(DataInput,Model)
            obj.DDMDataFit = DataInput;
            obj.FitResults = [];
            obj.Model = Model;
        end
        
        function DDMDataFitted = FitDDMData(obj)
            
            DDMFitOutput = table({},{},'VariableNames',{'Fit','CoeffVals'});

            for i=1:size(obj.DDMDataFit,1)
                time = obj.DDMDataFit.Time(i);
                DDMSignal = obj.DDMDataFit.DDMSignalValue(i);
                try
                f=  {fit( time{1,1}', DDMSignal{1,1}'./max( DDMSignal{1,1}),obj.Model)};
                DDMFitOutput.Fit(i) = f;
                DDMFitOutput.CoeffVals(i) = {coeffvalues(f{1,1})};
                catch
                DDMFitOutput.Fit(i) = {'Could Not Compute'};
                DDMFitOutput.CoeffVals(i) = {'Could Not Compute'};                   
                end
            end 
            DDMDataFitted = [obj.DDMDataFit, DDMFitOutput];
            obj.FitResults = DDMDataFitted;
        end
        
        function ShowFitResult(obj)
            colors = jet(size(obj.DDMDataFit,1));
            AllParametersInFunctionOFQ = [];
             for i=1:size(obj.DDMDataFit,1)
                try
                time =obj.FitResults.Time(i);
                DDMSignal = obj.FitResults.DDMSignalValue(i);
                f = obj.FitResults.Fit(i);
                CFVals = obj.FitResults.CoeffVals(i);
                QVector = obj.FitResults.QVector(i);
                
                time = time{1,1};
                DDMSignal = DDMSignal{1,1};
                f = f{1,1};
                CFVals = CFVals{1,1};
                if ~ischar(CFVals)
                     AllParametersInFunctionOFQ = [ AllParametersInFunctionOFQ ; [QVector , CFVals]];
                     NumPlots =numel(CFVals);
                end
                subplot(2,NumPlots,[1 NumPlots])
                xAxis = time(1):(time(end)/(5*length(time))):time(end);
                
                scatter(time,DDMSignal,50,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
                hold on
                plot(xAxis,f(xAxis).*max(DDMSignal),'Color',colors(i,:));
                catch
                     
                end

             end 
             
            set(gca,'LineWidth',2)
            set(gca,'FontSize',24)
            box on
            set(gca,'YScale','log')
            for i=1:NumPlots
                subplot(2,NumPlots, NumPlots+i)
                hold on
                box on
                scatter(AllParametersInFunctionOFQ(:,1),AllParametersInFunctionOFQ(:,i+1),50,'MarkerFaceColor','k','MarkerEdgeColor','k')
                set(gca,'LineWidth',2)
                set(gca,'FontSize',24)
                xlabel('Q-Vector')
                ylabel(['Parameter ' 'c_{' num2str(i) '}'])
            end
        end
        
        function ShowParameterResult(obj)
            
            
            
        end
    
end  
end