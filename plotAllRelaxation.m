clear;
clc;
close all;

folder2Data = 'D:\Boris_Hongbo\2_after adding cytoD_5 time points';
conditions = {'Cell1' 'Cell2' 'Cell3' 'Cell4' 'Cell5'};
conditionsVar = {'A','B','C', 'D', 'E'};

xAxis = categorical(conditions);
xTime = [0 20 40 60 80];

%% Load

files = dir(folder2Data);
counter = 1;
data = struct();

c1 = 1;
c2 = 1;
c3 = 1;
c4 = 1;
c5 = 1;
for i = 1:size(files,1)
    currentRes = [files(i).folder filesep files(i).name filesep 'results.mat'];
    if contains(files(i).name,'Cell') && files(i).isdir && isfile(currentRes)
        tmp = load([files(i).folder filesep files(i).name filesep 'results.mat']);
        
        
        if contains(files(i).name,conditions{1})
            data(c1).C1 = tmp.res;
            c1 = c1+1;     
        end
        
        if contains(files(i).name,conditions{2})
            data(c2).C2 = tmp.res;
            c2 = c2+1;            
        end
        
        if contains(files(i).name,conditions{3})
            data(c3).C3 = tmp.res;
            c3 = c3+1;            
        end
        
        if contains(files(i).name,conditions{4})
            data(c4).C4 = tmp.res;
            c4 = c4+1;            
        end
        
        if contains(files(i).name,conditions{5})
            data(c5).C5 = tmp.res;
            c5 = c5+1;            
        end
        
    end
    
end

%% Intensity comparison
fieldN = fieldnames(data);
nData = numel(fieldN);
nCond = length(xTime);
cellInt = zeros(nData,nCond);
polInt  = zeros(nData,nCond);
for i = 1:numel(fieldN)
    currF = ['C' num2str(i)];
    cellInt(i,:) = data.(currF).stats.cellInt(1:nCond);
    polInt(i,:)  = data.(currF).stats.polInt(1:nCond);
    
end
% shadow error bar
figure
subplot(1,2,1)
Plotting.shadedErrorBar(xTime,median(cellInt,1),std(cellInt,1))
title('Cell intensity')
xlim([xTime(1) xTime(end)]);
axis square
box on
subplot(1,2,2)
hold on
Plotting.shadedErrorBar(xTime,median(polInt,1),std(polInt,1))
title('Polymer intensity')
axis square
box on
xlim([xTime(1) xTime(end)]);

% Decay plot

figure
for i = 1:nCond
subplot(1,2,1)
    hold on
    plot(xTime,cellInt(i,:)/max(cellInt(i,:)))
    title('Cell intensity')
    xlim([xTime(1) xTime(end)]);
    axis square
    box on
    subplot(1,2,2)
    hold on
    plot(xTime,polInt(i,:)/max(polInt(i,:)))
    title('Polymer intensity')
    axis square
    box on
    xlim([xTime(1) xTime(end)]);
end

%% Volume comparison
fieldN = fieldnames(data);
nData = numel(fieldN);
nCond = length(xTime);
cellVol = zeros(nData,nCond);
polVol  = zeros(nData,nCond);
for i = 1:numel(fieldN)
    currF = ['C' num2str(i)];
    cellVol(i,:) = data.(currF).stats.cellVol(1:nCond);
    polVol(i,:)  = data.(currF).stats.polVol(1:nCond);
    
end
% shadow error bar
figure
subplot(1,2,1)
Plotting.shadedErrorBar(xTime,median(cellVol,1),std(cellVol,1))
title('Cell Volume')
xlim([xTime(1) xTime(end)]);
axis square
box on
subplot(1,2,2)
hold on
Plotting.shadedErrorBar(xTime,median(polVol,1),std(polVol,1))
title('Polymer Volume')
axis square
box on
xlim([xTime(1) xTime(end)]);

% Decay plot

figure
for i = 1:nCond
subplot(1,2,1)
    hold on
    plot(xTime,cellVol(i,:)/max(cellVol(i,:)))
    title('Cell Volume')
    xlim([xTime(1) xTime(end)]);
    axis square
    box on
    subplot(1,2,2)
    hold on
    plot(xTime,polVol(i,:)/max(polVol(i,:)))
    title('Polymer Volume')
    axis square
    box on
    xlim([xTime(1) xTime(end)]);
end

