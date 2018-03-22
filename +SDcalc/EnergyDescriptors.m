function [ vals, names ] = EnergyDescriptors( contour )
%NEWDESCRIPTOR Summary of this function goes here
%   Detailed explanation goes here

names = {'Elastic_Energy'; 'Bending_Energy'};
vals = cell(length(names),1);

if isempty(contour)
    return
else
    assert(size(contour,1) == 2, 'contour must be a [2xn] matrix' );
    assert(size(contour,2) > 2, 'a contour must constist of more than 2 points' );
    x = contour (1,:);
    y = contour (2,:);
    x = x - mean(x);
    y = y - mean(y);
    Xd1 = diff(x,1);
    Yd1 = diff(y,1);
    
    Xd2 = diff(x,2);
    Yd2 = diff(y,2);
    
    elas_ener = 0.5*(sum(Xd1.^2) + sum(Yd1.^2));
    bend_ener = 0.5*(sum(Xd2.^2) + sum(Yd2.^2));

% % %   % sanity check
% % %     ang = 30; % deg
% % %     ang = ang * pi / 180;
% % %     
% % %     R = [cos(ang), -sin(ang); sin(ang), cos(ang)];
% % %     
% % %     B2 = (R * B')';
% % %     
% % %     Xd1 = diff(B2(:,1),1);
% % %     Yd1 = diff(B2(:,2),1);
% % %     
% % %     Xd2 = diff(B2(:,1),2);
% % %     Yd2 = diff(B2(:,2),2);
% % %     
% % %     elas_ener2 = 0.5*(sum(Xd1.^2) + sum(Yd1.^2));
% % %     bend_ener2 = 0.5*(sum(Xd2.^2) + sum(Yd2.^2));
    
    v = [elas_ener; bend_ener];
    vals = num2cell(v);

end

end

