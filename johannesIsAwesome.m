clear
close all
clc

load('N:\Project Z-stacks\1K- 0.5 mg per ml\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps.mat')
porescold = allData;
load('Z:\Leica TCS SP8X\Leica PC\Johannes\2018\27052018\tif\sted\poreProps.mat')
poresheated = allData;
clear allData;
%% Plotting
figure(1)
totalV = porescold.Area;
totalV = totalV(:);
[CDF,CCDF] = Plotting.getCDF(totalV);
pcold.CCDF = CCDF;

totalV = poresheated.Area;
totalV = totalV(:);
[CDF,CCDF] = Plotting.getCDF(totalV);
pheated.CCDF = CCDF;

clear CDF CCDF totalV
figure()
CCDF = pcold.CCDF;
plot(CCDF.x,CCDF.y,'g')
hold on
CCDF = pheated.CCDF;
plot(CCDF.x,CCDF.y,'b')
hold off


a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF-plot for area of pores','confocal vs sted'})
xlabel('Pore size (micrometer^2)')
ylabel('CCDF - prob [0-1]')
a.FontSize = 14;
ylim([1e-6 1])
xlim([0.1 1e5])
legend({'confocal','sted'})

%%
figure(2)
subplot(2,1,1)
histogram(porescold.Solidity,0.3:0.01:0.99,'normalization','probability')
title({'Solidity of the pores',' PIC hydrogel (1K; 0.75 mg/ml) '})
xlabel('Solidity')
ylabel('normalised probability')
subplot(2,1,2)
histogram(poresheated.Solidity,0.3:0.01:0.99,'normalization','probability')
% subplot(3,1,3)
% histogram(porescold.Solidity,0.3:0.01:0.99,'normalization','probability')

%%
figure(3)
subplot(2,1,1)
histogram(porescold.Eccentricity,0.1:0.01:0.99,'normalization','probability')
subplot(2,1,2)
histogram(poresheated.Eccentricity,0.1:0.01:0.99,'normalization','probability')
% subplot(3,1,3)
% histogram(pores10.Eccentricity,0.1:0.01:0.99,'normalization','probability')