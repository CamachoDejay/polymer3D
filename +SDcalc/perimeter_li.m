function [vals, names] = perimeter_li( contour )
%x and y must be non-scalar vectors of equal sizes defining a polygonal
%region in clockwise or counterclockwise sequence. 
%   Detailed explanation goes here

names = {'Perimeter'};
vals = cell(length(names),1);

if isempty(contour)
    return
else
assert(size(contour,1) == 2, 'contour must be a [2xn] matrix' );
assert(size(contour,2) > 2, 'a contour must constist of more than 2 points' );

x = contour (1,:);
y = contour (2,:);
[m1,n1] = size(x); [m2,n2] = size(y);
n = max(m1,n1);
if logical (sum([m1,n1] ~= [m2,n2])) || min(m1,n1) ~= 1 || n <= 1
 error('Args must be equal-sized non-scalar vectors')
end
x = x(:); y = y(:);

if contour (:,1) == contour (:,end)
   % then contour is closed xn = x1
   x(end) = [];
   y(end) = [];
   n = n-1;
   % now I contour is open
end

x2 = [x(2:n);x(1)];
y2 = [y(2:n);y(1)];

x_diff  = x - x2;
y_diff  = y - y2;


eu_dis  = (x_diff.^2 + y_diff.^2).^0.5;

perim = sum(eu_dis(:));

vals(1) = {perim};


end
end



