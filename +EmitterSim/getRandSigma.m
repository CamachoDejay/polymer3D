function [ sigX, sigY ] = getRandSigma( setupPSFWidth, maxElip, emitter_number )
%GETRANDSigma random width of emitters 

%   Detailed explanation goes here
sigX = zeros(emitter_number,1);
sigY = zeros(emitter_number,1);

switch emitter_number
    case 1
        test = rand(1);
        if test>0.5
            sigX(1:round(emitter_number/2)) = setupPSFWidth+rand(1,round(emitter_number/2))*setupPSFWidth*(maxElip-1);
            sigY(1:round(emitter_number/2)) = setupPSFWidth+rand(1,round(emitter_number/2))*0.1;
        else
            sigX(1:round(emitter_number/2)) = setupPSFWidth+rand(1,round(emitter_number/2))*0.1;
            sigY(1:round(emitter_number/2)) = setupPSFWidth+rand(1,round(emitter_number/2))*setupPSFWidth*(maxElip-1);
        end
        
    otherwise
        sigX(1:round(emitter_number/2)) = setupPSFWidth+rand(1,round(emitter_number/2))*setupPSFWidth*(maxElip-1);
        sigY(1:round(emitter_number/2)) = setupPSFWidth+rand(1,round(emitter_number/2))*0.1;

        sigX(round(emitter_number/2)+1:end) = setupPSFWidth+rand(1,round(emitter_number/2))*0.1;
        sigY(round(emitter_number/2)+1:end) = setupPSFWidth+rand(1,round(emitter_number/2))*setupPSFWidth*(maxElip-1);
end
end

