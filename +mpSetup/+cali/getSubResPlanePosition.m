function [ zFocus, Fit ] = getSubResPlanePosition(focusMet, ZPos)

guess.mu = sum(ZPos.*focusMet)./sum(focusMet);
guess.sig = 2;
zFocus = zeros(1,size(focusMet,2));
for k=1:size(focusMet,2)
    [out,Fit] = SimpleFitting.gauss1D(focusMet(:,k), ZPos,guess);
    zFocus(k) = out(2);
end
