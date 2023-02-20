function [ zFocus, Fit ] = getSubResPlanePosition(focusMet, ZPos)

[val,id] = max(focusMet);
guess.mu = ZPos(id);
guess.sig = 1;
zFocus = zeros(1,size(focusMet,2));

% Find the half max value.
halfMax = (min(focusMet) + max(focusMet)) / 2;


% Find where the data first drops below half the max.
beforeMax = flipud(focusMet(1:id));
afterMax  = focusMet(id:end);

id1 = find(beforeMax<=halfMax,1,'first');
index1 = id-id1;
%index1 = find(focusMet >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
id2 = find(afterMax<=halfMax,1,"first");
index2 = id+id2;

%index2 = find(focusMet>= halfMax, 1, 'last');

fwhmx = ZPos(index2) - ZPos(index1);

delta = fwhmx/2;
range = and(ZPos>ZPos(id)-delta,ZPos<ZPos(id)+delta);

ZPos2Use = ZPos(range);
focMet2Use = focusMet(range,:);

[out,Fit] = SimpleFitting.gauss1D(focMet2Use, ZPos2Use,guess);
zFocus = out(2);
Fit(:,1) = ZPos2Use;
Fit(:,2) = SimpleFitting.gaussian(out,ZPos2Use);

