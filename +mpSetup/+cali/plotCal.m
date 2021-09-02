function plotCal(ch1,ch2,ROI)

figure
for i= 1:4
subplot(2,4,i)
imagesc(ch1(ROI(i,2):ROI(i,2)+ROI(i,4)-1,ROI(i,1):ROI(i,1)+ROI(i,3)-1));


end

for i= 1:4
subplot(2,4,4+i)
imagesc(ch2(ROI(4+i,2):ROI(4+i,2)+ROI(4+i,4)-1,ROI(4+i,1):ROI(4+i,1)+ROI(4+i,3)-1));


end
end