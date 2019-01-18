clear
close all
clc

load('N:\Project Z-stacks\New stacks\1K_0.5mg_30 nm_sample 1\Tif_1K_0,5mg_30 nm_sample 1\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps.mat')
porescold = allDataAdapt;
load('N:\Project Z-stacks\New stacks\5k_0.5 mg 30 nm_Sample 1\tif_5k_0.5 mg 30 nm_Sample 1\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps.mat')
poresheated = allDataAdapt;
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