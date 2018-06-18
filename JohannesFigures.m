% The purpose of this program is to: 
clear
close all
clc

%% Load block
% Load first file
tmp = load('N:\Project pnipam\tif cold 1\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps.mat');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
Table05 = tmp.(names{1});
clear tmp names
% Load second file
tmp = load('N:\Project pnipam\tif hot 1\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps.mat');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
Table075 = tmp.(names{1});
clear tmp names
% Loqd third file
tmp = load('N:\Project pnipam\tif cold 2\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps.mat');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
Table1 = tmp.(names{1});
clear tmp names

%% Output block


% ask for area an put that into totalV
totalV = Table05.Area;
% make sure is a col vector
totalV = totalV(:);
% calculate CDF
[CDF,CCDF] = Plotting.getCDF(totalV);
% store CCDF into p05
p05_CCDF = CCDF;

% qsk for area an put that into totalV
totalV = Table075.Area;
% make sure is a col vector
totalV = totalV(:);
% calculate CDF
[CDF,CCDF] = Plotting.getCDF(totalV);
% store CCDF into p05
p075_CCDF = CCDF;

% qsk for area an put that into totalV
totalV = Table1.Area;
% make sure is a col vector
totalV = totalV(:);
% calculate CDF
[CDF,CCDF] = Plotting.getCDF(totalV);
% store CCDF into p05
p1_CCDF = CCDF;
clear CDF CCDF totalV

%%
figure()

CCDF = p05_CCDF;
plot(CCDF.x,CCDF.y,'g')
hold on

CCDF = p075_CCDF;
plot(CCDF.x,CCDF.y,'b')

CCDF = p1_CCDF;
plot(CCDF.x,CCDF.y,'r')
hold off

a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF-plot for area of pores','PIC/pnipam; hydrogel 0.5 mg/ml heat response'})
xlabel('Pore size (micrometer^2)')
ylabel('CCDF - prob [0-1]')
a.FontSize = 14;
ylim([1e-6 1])
xlim([0.1 1e5])
legend({'cold 1','hot 1','cold 2'})

%%
figure(2)
subplot(3,1,1)
histogram(pores05.Solidity,0.3:0.01:0.99,'normalization','probability')
title({'Solidity of the pores',' PIC hydrogel (1K; 0.75 mg/ml) '})
xlabel('Solidity')
ylabel('normalised probability')
subplot(3,1,2)
histogram(pores07.Solidity,0.3:0.01:0.99,'normalization','probability')
subplot(3,1,3)
histogram(pores10.Solidity,0.3:0.01:0.99,'normalization','probability')


figure(3)
subplot(3,1,1)
histogram(pores05.Eccentricity,0.1:0.01:0.99,'normalization','probability')
subplot(3,1,2)
histogram(pores07.Eccentricity,0.1:0.01:0.99,'normalization','probability')
subplot(3,1,3)
histogram(pores10.Eccentricity,0.1:0.01:0.99,'normalization','probability')