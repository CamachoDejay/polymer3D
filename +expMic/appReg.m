function [contour] = appReg(contour,angle,scaling)
%APPREG applies simple registration to contour, only rotation and scaling
%are considered, angular rotation must be in degrees.
%   Detailed explanation goes here

    assert(size(contour,1) == 2, 'contour must be a [2xn] matrix' );
    assert(size(contour,2) > 2, 'contour must contain more than 2 points');
    % change to radians
    T = angle*pi/180;
    % rotation matrix
    rotMat = [cos(T), -sin(T); sin(T), cos(T)];
    % so the matrix operation works
    contour = contour';
    % rotation
    contour = contour*rotMat;
    % scaling
    contour = contour.*scaling;
    % back to the old 2xn matrix
    contour = contour';
end

