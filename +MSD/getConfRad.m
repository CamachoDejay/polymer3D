% Obtain the confinement volume from msd trace. By fitting an exponential
% trace A0-A0*exp(X) we can calculate the square root of A0 (tranform from
% um^2 to um = radius, in which the trace is confined. 
% 
% WARNING: There is no check that the motion is actually confined motion,
% usually one wants to determine alpha and ensure that it is significantly
% lower than 1 before assuming that the motion is confined. This code will
% always output a value, so one need to be mindful before using it.
% The data and the fit are plotted so it is easy to keep a critical eye on
% the fit
%
% INPUT:
%   -  msd: msd curve of the data presented as a single column vector
%   -  fitRange: [0 1]=> fraction of the data to be used for fitting
%      starting from beginning
% 
% OUTPUT:
%   - rConf: Confinement radius = sqrt(y0+A0), unit is the same as the
%            input unit of msd;

function [rConf] = getConfRad(msd,fitRange,expTime)
     
    domain = (1:round(fitRange*length(msd)))*expTime;
    toFit  = msd(1:round(fitRange*length(msd)));
    initGuess = [max(toFit)-min(toFit),min(toFit),mean(toFit)];
    
    [fPar,resnorm,res] = lsqcurvefit(@invExp,initGuess,...
            domain(:),toFit(:));
    
     fit = invExp(fPar,domain);
     
     figure
     plot(domain,toFit)
     hold on
     plot(domain,fit)
     
     rConf = sqrt(fPar(1)+fPar(2));

end


% Function used for fitting A-A*exp(-X)+y0;
function F = invExp(x,data)
    
    F = x(1) * ( 1 - exp( - x(3)*data/x(1)))+x(2);


end