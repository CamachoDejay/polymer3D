        function [binnedData] = qBinning(data2Bin,nPoint)
            
            %small function that bins the data using quantiles
            nQuant = linspace(0,1,length(data2Bin)/nPoint);
            binnedData = quantile(data2Bin,nQuant);
            
        end