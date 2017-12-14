function [ blank_bg ] = gauss_noise( dim, mean_val, FWHM, d_range )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% mean_val = 1000;
% FWHM = 100;
% d_range   = 65000;
if strcmp(d_range,'uint16')
    d_range_val = 65535;
else
    error('Problems with dynamic range')
end

M  = mean_val/d_range_val;
FW = FWHM/d_range_val;
S  = FW / (2*((2*log(2))^0.5)); % std
V = S^2;

blank = zeros(dim);
blank_bg = imnoise(blank,'gaussian',M,V);
blank_bg = blank_bg.*d_range_val;

if strcmp(d_range,'uint16')
    blank_bg = uint16(blank_bg);
end

end

