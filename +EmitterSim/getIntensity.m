function [ em_int ] = getIntensity( model, em_n )
%GETINTENSITY calculates the intensity of each emitter in the image using
%the desired model
%   Detailed explanation goes here
n = model.name;

if strcmp(n,'normal')
    
    em_mean_int = model.mean_int;
    em_int_sigma = model.sigma;

    em_int = randn(em_n,1);
    em_int = repmat(em_mean_int,size(em_int)) + em_int.*em_int_sigma;    

end
end

