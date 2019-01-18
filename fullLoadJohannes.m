clear
close all
clc

%Explanation:
%This code aim to be used to load data all data so one can use
%fullPlotJohannes to plot the loaded data.

%There is one executable block per Sample In which Log - Histogram and CCDF are
%calculated keeping track of individual stack and having all stack
%together.

%% User Input
stack2Use = 5; % Take the last 5 stacks
numBin = 20; %Number of Bins
% Load block 1K
%% S1_1K
tmp = load('N:\Project Z-stacks\New stacks\1K_0.25mg_30 nm_sample 1\tif_1K_0.25 mg-30 nm_Sample 1\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S1_1K.Table025 = tmp.(names{1});
clear tmp names
% Load second file
tmp = load('N:\Project Z-stacks\New stacks\1K_0.5mg_30 nm_sample 1\Tif_1K_0,5mg_30 nm_sample 1\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S1_1K.Table050 = tmp.(names{1});
clear tmp names
% Loqd third file
tmp = load('N:\Project Z-stacks\New stacks\1K_1 mg_30 nm_sample 1\TIF_1K_1 mg_0_30 nm_Sample1\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S1_1K.Table100 = tmp.(names{1});
clear tmp names

%CDF
name = fieldnames(S1_1K);
numField = length(name);
for i = 1 : numField 
    
    S1_1K.(name{i})(1).CCDF = [];
    allArea =[];
    
    for j = 1: stack2Use
        
        area = S1_1K.(name{i})(end+1-j).Data.Area;
        area = area(:);
        [CDF,CCDF] = Plotting.getCDF(area);
        CCDF = [CCDF.x CCDF.y];
        S1_1K.(name{i})(end+1-j).CCDF = CCDF;
        allArea = [allArea ; area];
        
    end
    
    S1_1K.(name{i})(end+1-j).allArea = allArea;
    [CDF,CCDF] = Plotting.getCDF(allArea(:));
    CCDF = [CCDF.x CCDF.y];
    S1_1K.(name{i})(end+1-j).allCCDF = CCDF;
    
    %for Histogram
    [midBin, Freq, LFreq] = Plotting.lnbin(allArea,numBin);
    
    S1_1K.(name{i})(end+1-j).midBin = midBin;
    S1_1K.(name{i})(end+1-j).Freq = Freq;
    S1_1K.(name{i})(end+1-j).LFreq = LFreq;
    
    if length(S1_1K.(name{i}))-j>0
        
        S1_1K.(name{i})(1:length(S1_1K.(name{i}))-j) = [];
    
    end
    
end

%% S2_1K
tmp = load('N:\Project Z-stacks\New stacks\1K_0,25mg_30 nm_sample 2\TIF_1K_0.25 mg_30 nm_sample2\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S2_1K.Table025 = tmp.(names{1});
clear tmp names
% Load second file
tmp = load('N:\Project Z-stacks\New stacks\1K_0,5mg_30 nm_sample 2\tif_1k_0.5 mg_30 nm_Sample 2\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S2_1K.Table050 = tmp.(names{1});
clear tmp names
% Loqd third file
tmp = load('N:\Project Z-stacks\New stacks\1K_1 mg_30 nm_sample 2\tif_1K_1 mg_30 nm_Sample 2\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S2_1K.Table100 = tmp.(names{1});
clear tmp names

%CDF
name = fieldnames(S2_1K);
 
for i = 1 : numField 
    
    S2_1K.(name{i})(1).CCDF = [];
    allArea =[];
    
    for j = 1: stack2Use
        
        area = S2_1K.(name{i})(end+1-j).Data.Area;
        area = area(:);
        [CDF,CCDF] = Plotting.getCDF(area);
        CCDF = [CCDF.x CCDF.y];
        S2_1K.(name{i})(end+1-j).CCDF = CCDF;
        allArea = [allArea ; area];
        
    end
    
    S2_1K.(name{i})(end+1-j).allArea = allArea;
    [CDF,CCDF] = Plotting.getCDF(allArea(:));
    CCDF = [CCDF.x CCDF.y];
    S2_1K.(name{i})(end+1-j).allCCDF = CCDF;
    
    %for Histogram
    [midBin, Freq, LFreq] = Plotting.lnbin(allArea,numBin);
    
    S2_1K.(name{i})(end+1-j).midBin = midBin;
    S2_1K.(name{i})(end+1-j).Freq = Freq;
    S2_1K.(name{i})(end+1-j).LFreq = LFreq;
    
    if length(S2_1K.(name{i}))-j>0
        
        S2_1K.(name{i})(1:length(S2_1K.(name{i}))-j) = [];
    
    end
end

%Load block 5K

%% S1_5K
tmp = load('N:\Project Z-stacks\New stacks\5K_0.25 mg_30 nm_Sample 1\tif_5K_0.25 mg_30 nm_Sample 1\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S1_5K.Table025 = tmp.(names{1});
clear tmp names
% Load second file
tmp = load('N:\Project Z-stacks\New stacks\5k_0.5 mg 30 nm_Sample 1\tif_5k_0.5 mg 30 nm_Sample 1\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S1_5K.Table050 = tmp.(names{1});
clear tmp names
% Loqd third file
tmp = load('N:\Project Z-stacks\New stacks\5K_1 mg_30nm_Sample 1\tif_5K_1 mg 30 nm_Sample 1\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S1_5K.Table100 = tmp.(names{1});
clear tmp names

%CDF
name = fieldnames(S1_5K);
 
for i = 1 : numField 
    
    S1_5K.(name{i})(1).CCDF = [];
    allArea =[];
    
    for j = 1: stack2Use
        
        area = S1_5K.(name{i})(end+1-j).Data.Area;
        area = area(:);
        [CDF,CCDF] = Plotting.getCDF(area);
        CCDF = [CCDF.x CCDF.y];
        S1_5K.(name{i})(end+1-j).CCDF = CCDF;
        allArea = [allArea ; area];
        
    end
    
    S1_5K.(name{i})(end+1-j).allArea = allArea;
    [CDF,CCDF] = Plotting.getCDF(allArea(:));
    CCDF = [CCDF.x CCDF.y];
    S1_5K.(name{i})(end+1-j).allCCDF = CCDF;
    
    %for Histogram
    [midBin, Freq, LFreq] = Plotting.lnbin(allArea,numBin);
    
    S1_5K.(name{i})(end+1-j).midBin = midBin;
    S1_5K.(name{i})(end+1-j).Freq = Freq;
    S1_5K.(name{i})(end+1-j).LFreq = LFreq;
    
    if length(S1_5K.(name{i}))-j>0
        
        S1_5K.(name{i})(1:length(S1_5K.(name{i}))-j) = [];
    
    end
    
end

%% S2_5K

tmp = load('N:\Project Z-stacks\New stacks\5K_0.25mg_30 nm_Sample 2\tif_5K_0.25mg_30 nm_Sample 2\correct stack\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S2_5K.Table025 = tmp.(names{1});
clear tmp names
% Load second file
tmp = load('N:\Project Z-stacks\New stacks\5k_0.5 mg 30 nm_Sample 2\tif_5K_0.5mg_30 nm_Sample 2\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S2_5K.Table050 = tmp.(names{1});
clear tmp names
% Loqd third file
tmp = load('N:\Project Z-stacks\New stacks\5K_1 mg_30nm_Sample 2\tif_5K_1 mg_30nm_Sample 2\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
S2_5K.Table100 = tmp.(names{1});
clear tmp names

%CDF
name = fieldnames(S2_5K);
 
for i = 1 : numField 
    
    S2_5K.(name{i})(1).CCDF = [];
    allArea =[];
    
    for j = 1: stack2Use
        
        area = S2_5K.(name{i})(end+1-j).Data.Area;
        area = area(:);
        [CDF,CCDF] = Plotting.getCDF(area);
        CCDF = [CCDF.x CCDF.y];
        S2_5K.(name{i})(end+1-j).CCDF = CCDF;
        allArea = [allArea ; area];
        
    end
    
    S2_5K.(name{i})(end+1-j).allArea = allArea;
    [CDF,CCDF] = Plotting.getCDF(allArea(:));
    CCDF = [CCDF.x CCDF.y];
    S2_5K.(name{i})(end+1-j).allCCDF = CCDF;
    
    %for Histogram
    [midBin, Freq, LFreq] = Plotting.lnbin(allArea,numBin);
    
    S2_5K.(name{i})(end+1-j).midBin = midBin;
    S2_5K.(name{i})(end+1-j).Freq = Freq;
    S2_5K.(name{i})(end+1-j).LFreq = LFreq;
    
    if length(S2_5K.(name{i}))-j>0
        
        S2_5K.(name{i})(1:length(S2_5K.(name{i}))-j) = [];
    
    end
    
end
disp('============> DONE <============');