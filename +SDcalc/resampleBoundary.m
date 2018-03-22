function [ boundary_RS ] = resampleBoundary( boundary, number_of_points )
% Resamples a boundary using smoothing spline. Boundary is expected to ahve
% the same format as that given by the function bwboundaries, where the
% first and last element are the same. 
%   Detailed explanation goes here

% check that boundaries come with expected format
assert(size(boundary,2) == 2, 'Boundary must consist of a 2 col matrix')
assert(size(boundary,1) >= 3, 'Boundary must consist of more than 2 points')
if boundary(1,:) == boundary(end,:)
    b = boundary(1:end-1,:);
%     disp('Working on boundary')
else
    warning('Boundary not given in expected format')
    boundary_RS = [];
    return
end

% initial number of unique points in boundary
points_n_i  = length(b);

% find how much we have to extend the data to avoid edge efects
if number_of_points > points_n_i;
    margin = round(points_n_i/10);
else
    margin = round(number_of_points/10);
end

% we extend the contour using known data, this is posible because we are
% working with the cotour of a closed shape
extend_1   =  b(end-margin+1:end,:);
extend_2   =  b(1:margin,:);
extended_b = [extend_1; b; extend_2];

% preparing data for fitting
x_vals     = extended_b(:,1);
y_vals     = extended_b(:,2);
points_n_e = length(x_vals);
c          = (1:points_n_e)';

% fitting
fit_type    = 'smoothingspline';%'linearinterp';%
fitobject_x   = fit(c,x_vals,fit_type);
fitobject_y   = fit(c,y_vals,fit_type);

% obtaining fitted values
n_points      = number_of_points +1 ;
c_i           = c(margin+1);
c_f           = c(points_n_i+margin+1);
c_fit         = linspace(c_i,c_f,n_points); c_fit = c_fit(:);
x_vals_fit    = feval(fitobject_x,c_fit);
y_vals_fit    = feval(fitobject_y,c_fit);

% getting the right range of values and closing the boundary to keep
% formatting agreement
tmp  = [x_vals_fit y_vals_fit];
boundary_RS = tmp(1:end-1,:);
boundary_RS = [boundary_RS; boundary_RS(1,:)];

end

