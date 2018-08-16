function [ zFocus, Fit ] = getSubResPlanePosition(focusMet, ZPos)

zFocus = zeros(1,size(focusMet,2));
for k=1:size(focusMet,2)
    [out,Fit] = SimpleFitting.gauss1D(focusMet(:,k), ZPos);
    zFocus(k) = out(2);
end
