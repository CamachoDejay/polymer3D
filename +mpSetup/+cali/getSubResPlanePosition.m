function [ zFocus, Fit ] = getSubResPlanePosition(focusMet, ZPos)

[val,id] = max(focusMet);
guess.mu = ZPos(id);
guess.sig = 1;
zFocus = zeros(1,size(focusMet,2));
delta = 1.5;
range = and(ZPos>ZPos(id)-delta,ZPos<ZPos(id)+delta);

ZPos2Use = ZPos(range);
focMet2Use = focusMet(range,:);

[out,Fit] = SimpleFitting.gauss1D(focMet2Use, ZPos2Use,guess);
zFocus = out(2);
Fit(:,1) = ZPos2Use;
Fit(:,2) = SimpleFitting.gaussian(out,ZPos2Use);


