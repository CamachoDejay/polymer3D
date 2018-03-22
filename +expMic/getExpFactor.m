function [linExp, reg] = getExpFactor(Pre,Post)
%GETEXPFACTOR get linear expansion factor based on the contours of the Pre
%and Post expanded structures.
%   Detailed explanation goes here
assert(isstruct(Pre), 'Bad input, Pre must be the right structure')
assert(isstruct(Post), 'Bad input, Pre must be the right structure')

fn = fieldnames(Pre);
expectedFn = {'rec';'contour';'centroid';'area';'boundingBox';'length';'ang';'contNorm'};
assert(all(strcmp(fn,expectedFn)), 'Pre field names are not ok')

fn = fieldnames(Post);
expectedFn = {'rec';'contour';'centroid';'area';'boundingBox';'length';'ang';'contNorm'};
assert(all(strcmp(fn,expectedFn)), 'Post field names are not ok')

% linear scaling factor
scaling = Pre.length/Post.length;

% orienttion of the structures
ang1 = Pre.ang;
ang2 = Post.ang;

% we have to check which one is the right angle
theta = [ang2-ang1, ang2-ang1+180];
% mean deviation storage
meanD = nan(2,1);
% registered contour storage
regCont = zeros([size(Post.contNorm), 2]);
% Fixted contour
fixed = Pre.contNorm;
for i = 1:2
    % apply registration
    [regPost] = expMic.appReg(Post.contNorm,theta(i),scaling);
    % registered contour
    regCont(:,:,i) = regPost;
    % mean deviation based on closest neighbour
    [~, D] = knnsearch(fixed',regPost');
    meanD(i) = mean(D);

end

% find the best fit
[errVal,idx] = min(meanD);
% store best values
theta = theta(idx);
regCont = regCont(:,:,idx);

% storing registration information
reg.scaling = scaling;
reg.rot = theta;
reg.contour = regCont;
reg.error = errVal;

% output for linear expansion factor
linExp = 1/scaling;
end

