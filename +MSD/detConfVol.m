function [rConf] = detConfVol(msd,expTime)
     
    domain = (1:length(msd))*expTime;
    initGuess = [median(msd),max(msd),mean(msd)];
    [fPar,resnorm,res] = lsqcurvefit(@invExp,initGuess,...
            domain(:),msd(:));
    
     fit = invExp(fPar,domain);
     
     figure
     plot(domain,msd)
     hold on
     plot(domain,fit)
     
     rConf = fPar(1);

end





function F = invExp(x,data)
    
    F = x(1) * ( 1 - x(2)*exp( - 4*x(3)*data/x(1)));


end