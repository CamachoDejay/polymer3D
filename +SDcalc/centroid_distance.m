function [ c_distance ] = centroid_distance( contour )
%The centroid distance function -r- is expresed by the distance of the
%boundary points from the centroid (x_c,y_c) of the shape. r(t) where t is
%defined 0<t<T, T is the period of t. r(t) is invariant to translation. 
%   Detailed explanation goes here

if isempty(contour)
    % output generation
    c_distance.raw.r = [];
    c_distance.raw.t = [];
    c_distance.raw.theta = [];
    c_distance.fit.r = [];
    c_distance.fit.t = [];
    
else
    assert(size(contour,1) == 2, 'contour must be a [2xn] matrix' );
    assert(size(contour,2) > 2, 'a contour must constist of more than 2 points' );

    x = contour (1,:);
    y = contour (2,:);
    x = x(:); y = y(:);
    n = length(x);

    [x_c,y_c,~] = centroid_by_area( contour );

    if contour (:,1) == contour (:,end)
       % then contour is closed xn = x1
       x(end) = [];
       y(end) = [];
       n = n-1;
       % now I contour is open
    end

    %%%%%%%%%     stimation of the perimeter    %%%%%%%%%%%%%%%%%%%%
    x2 = [x(2:n);x(1)];
    y2 = [y(2:n);y(1)];

    x_diff  = x - x2;
    y_diff  = y - y2;
    % eucledian distances
    d_l  = (x_diff.^2 + y_diff.^2).^0.5; d_l = d_l(:);
    L = sum(d_l(:));
    clear x2 y2 x_diff y_diff
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%% calculation of r(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_diff  = y - y_c;
    x_diff  = x - x_c;
    % eucledian distance to the centroid
    r = (x_diff.^2 + y_diff.^2).^0.5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% calculation of the frequency  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = zeros(size(r));
    for i = 2:length(x)
        t (i,1) = sum(d_l(1:i-1)) / L;
    end
    t = t.*2*pi; % this is to have the information on a frequency scale that 
                 % makes sence. 2pi refers to the complete arc. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% calculation of angular information %%%%%%%%%%%%%%%%%%%%%%%%
    % this is important if I wat to recreate the shape in a easy way 

    theta = atan (y_diff./x_diff);
    indx = x_diff > 0 & y_diff < 0;
    theta(indx) = theta(indx) + 2*pi;
    indx = x_diff < 0 ;
    theta(indx) = theta(indx) + pi;
    theta = theta.*-1;  %due to the rotation convention
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%% interpolation of curve to make frequency values invariant to the
    %%%%%% the initial contour

    fit_type    = 'linearinterp';%'smoothingspline';
    fitobject   = fit(t,r,fit_type);
    n_freq_points = 1000;
    t_fit         = linspace(0,(2*pi),n_freq_points+1); t_fit = t_fit(:);
    r_fit      = feval(fitobject,t_fit);

    r_fit(end) = [];
    t_fit(end)    = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % output generation
    c_distance.raw.r = r;
    c_distance.raw.t = t;
    c_distance.raw.theta = theta;
    c_distance.fit.r = r_fit;
    c_distance.fit.t = t_fit;
end

end

