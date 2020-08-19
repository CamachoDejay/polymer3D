function [RMSD, fitVal] = minGauss1D(x,y,params)
    A = params(3);
    sigma = params(1);
    mu    = params(2);
    y0    = params(4);
    nVal = length(x);
    fitVal   = y0 + (A/sqrt(2*pi*sigma^2)).*exp((-(x-mu).^2)./(2*sigma^2));

    dev    = fitVal-y;
    devSqr = (dev).^2;

    RMSD = sqrt(sum(devSqr)/nVal);
