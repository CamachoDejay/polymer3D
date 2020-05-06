%% Function to Calculate viscosity based on Stokes-Einstein equation
% INPUT : 
%   -  D: Diffusion of particle in um^2/s
%   -  r: Radius of the particle in um
%
% OUTPUT:
%   -  n: Viscosity in centipoise
function [n] = getViscosity(D,r,T)

    n = ((1.380649*10^(-23)*T)/(6*pi*r*10^(-6)*D*10^(-12)))*1000;%in cp


end