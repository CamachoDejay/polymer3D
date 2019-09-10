% This function aim at determining which type of diffusion we have by
% determining the anomaly diffusion coefficient alpha. According to the
% following source : https://tinevez.github.io/msdanalyzer/tutorial/MSDTuto_confined.html#4
% This can be done by fitting msd curve in a log-log scale and obtaining
% the slope as alpha.
%
% if alpha >1 : super-diffusive regime (e.g. active transport, drift)
% if alpha = 1: brownian motion
% if alpha < 1: sub-diffusive to confine
%
% INPUT:
%   - msd: msd curve for the data
%   - expTime: exposure time, time between two frame
% OUPUT:
%   - alpha: abnormal diffusion coefficient

function [alpha] = detDiffTypeAlpha(msd,expTime)

    t = (1:length(msd))*expTime;
    toFit = log(msd);
    fitPar = fit(log(t(:)),toFit(:),'a*x+b');
    alpha = fitPar.a;
   
end