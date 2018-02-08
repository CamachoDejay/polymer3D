function [RMSD, fitVal] = fun2minGauss(rawData,params,doPlot)

N = rawData.N;
binCenter  = rawData.binCenter;
nVal   = length(binCenter);


sigma1 = params(1);
mu1    = params(2);
A1     = params(3);
bkg    = params(4);

fitVal   = ((A1/sqrt(2*pi*sigma1^2)).*exp((-(binCenter-mu1).^2)./(2*sigma1^2)))+bkg;

dev    = fitVal-N;
devSqr = (dev).^2;

RMSD = sqrt(sum(devSqr)/nVal);


if doPlot
    figure(69)
    plot(N,'k')
    hold on
    plot(fitVal,'b')
    plot(dev,'r')
    hold off
    title(['RMSD: ' num2str(RMSD)])
    shg
    drawnow
    pause(0.1)
end
end