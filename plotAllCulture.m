clear;
clc;
close all;

folder2Data = 'D:\Boris_Hongbo\1_with culturing time_0-24h';
conditions = {'0h' '2h' '7h' '15h' '24h'};
conditionsVar = {'A','B','C', 'D', 'E'};

xAxis = categorical(conditions);
xTime = [0 2 7 15 24];
%% Load

files = dir(folder2Data);
counter = 1;
data = struct();

c0 = 1;
c1 = 1;
c2 = 1;
c3 = 1;
c4 = 1;
for i = 1:size(files,1)
    currentRes = [files(i).folder filesep files(i).name filesep 'results.mat'];
    if contains(files(i).name,'Cell') && files(i).isdir && isfile(currentRes)
        tmp = load([files(i).folder filesep files(i).name filesep 'results.mat']);
        
        
        if contains(files(i).name(end-1:end),conditions{1})
            data(c0).T0 = tmp.res;
            c0 = c0+1;     
        end
        
        if contains(files(i).name(end-1:end),conditions{2})
            data(c1).T1 = tmp.res;
            c1 = c1+1;            
        end
        
        if contains(files(i).name(end-1:end),conditions{3})
            data(c2).T2 = tmp.res;
            c2 = c2+1;            
        end
        
        if contains(files(i).name(end-2:end),conditions{4})
            data(c3).T3 = tmp.res;
            c3 = c3+1;            
        end
        
        if contains(files(i).name(end-2:end),conditions{5})
            data(c4).T4 = tmp.res;
            c4 = c4+1;            
        end
        
    end
    
end


%% Intensity comparison
fieldN = fieldnames(data);
cellInt = zeros(numel(fieldN),3);
polInt  = zeros(numel(fieldN),3);

for i = 1:numel(fieldN)
    currF = ['T' num2str(i-1)];
    
    tmpCellInt = zeros(size(data,2),1);
    tmpPolInt  = zeros(size(data,2),1);
    for j = 1:size(data,2)
       
        tmpCellInt(j) = data(j).(currF).stats.cellInt;
        tmpPolInt(j)  = data(j).(currF).stats.polInt;
        
    end
    
    cellInt(i,1) = median(tmpCellInt);
    cellInt(i,2) = std(tmpCellInt);
    cellInt(i,3) = std(tmpCellInt);
    
    polInt(i,1) = median(tmpPolInt);
    polInt(i,2) = std(tmpPolInt);
    polInt(i,3) = std(tmpPolInt);
    
end

xAxis = reordercats(xAxis,{'0h' '2h' '7h' '15h' '24h'});
% figure
% subplot(1,2,1)
% bar(xAxis,cellInt(:,1))
% hold on
% er = errorbar(xAxis,cellInt(:,1),cellInt(:,2),cellInt(:,3));
% er.Color =[0 0 0];
% er.LineStyle = 'none';
% title('Cell intensity')
% axis square
% 
% subplot(1,2,2)
% bar(xAxis,polInt(:,1))
% hold on
% er = errorbar(xAxis,polInt(:,1),polInt(:,2),polInt(:,3));
% er.Color =[0 0 0];
% er.LineStyle = 'none';
% title('Polymer intensity')
% axis square

figure
subplot(1,2,1)
Plotting.shadedErrorBar(xTime,cellInt(:,1),cellInt(:,2))
title('Cell intensity')
xlim([xTime(1) xTime(end)]);
axis square
box on
subplot(1,2,2)
hold on
Plotting.shadedErrorBar(xTime,polInt(:,1),polInt(:,2))
title('Polymer intensity')
axis square
box on
xlim([xTime(1) xTime(end)]);

%% Volume comparison
fieldN = fieldnames(data);
cellVol = zeros(numel(fieldN),3);
polVol  = zeros(numel(fieldN),3);

for i = 1:numel(fieldN)
    currF = ['T' num2str(i-1)];
    
    tmpCellVol = zeros(size(data,2),1);
    tmpPolVol  = zeros(size(data,2),1);
    for j = 1:size(data,2)
       
        tmpCellVol(j) = data(j).(currF).stats.cellVol;
        tmpPolVol(j)  = data(j).(currF).stats.polVol;
        
    end
    
    cellVol(i,1) = median(tmpCellVol);
    cellVol(i,2) = std(tmpCellVol);
    cellVol(i,3) = std(tmpCellVol);
    
    polVol(i,1) = median(tmpPolVol);
    polVol(i,2) = std(tmpPolVol);
    polVol(i,3) = std(tmpPolVol);
    
end

xAxis = reordercats(xAxis,{'0h' '2h' '7h' '15h' '24h'});
% figure
% subplot(1,2,1)
% bar(xAxis,cellVol(:,1))
% hold on
% er = errorbar(xAxis,cellVol(:,1),cellVol(:,2),cellVol(:,3));
% er.Color =[0 0 0];
% er.LineStyle = 'none';
% title('Cell Volensity')
% axis square
% 
% subplot(1,2,2)
% bar(xAxis,polVol(:,1))
% hold on
% er = errorbar(xAxis,polVol(:,1),polVol(:,2),polVol(:,3));
% er.Color =[0 0 0];
% er.LineStyle = 'none';
% title('Polymer Volensity')
% axis square

figure
subplot(1,2,1)
Plotting.shadedErrorBar(xTime,cellVol(:,1),cellVol(:,2))
title('Cell Volume')
xlim([xTime(1) xTime(end)]);
axis square
box on
subplot(1,2,2)
hold on
Plotting.shadedErrorBar(xTime,polVol(:,1),polVol(:,2))
title('Polymer Volume')
axis square
box on
xlim([xTime(1) xTime(end)]);
%% Volume normalized
fieldN = fieldnames(data);
cellVol = zeros(numel(fieldN),3);
polVol  = zeros(numel(fieldN),3);

for i = 1:numel(fieldN)
    currF = ['T' num2str(i-1)];
    
    tmpCellVol = zeros(size(data,2),1);
    tmpPolVol  = zeros(size(data,2),1);
    normPolVol = zeros(size(data,2),1);
    
    for j = 1:size(data,2)
       
        tmpCellVol(j) = data(j).(currF).stats.cellVol;
        tmpPolVol(j)  = data(j).(currF).stats.polVol;
        normPolVol(j) = tmpPolVol(j)/tmpCellVol(j);
        
    end
    
    polVol(i,1) = median(normPolVol);
    polVol(i,2) = std(normPolVol);
    polVol(i,3) = std(normPolVol);
    
end

xAxis = reordercats(xAxis,{'0h' '2h' '7h' '15h' '24h'});
% figure
% subplot(1,2,1)
% bar(xAxis,cellVol(:,1))
% hold on
% er = errorbar(xAxis,cellVol(:,1),cellVol(:,2),cellVol(:,3));
% er.Color =[0 0 0];
% er.LineStyle = 'none';
% title('Cell Volensity')
% axis square
% 
% subplot(1,2,2)
% bar(xAxis,polVol(:,1))
% hold on
% er = errorbar(xAxis,polVol(:,1),polVol(:,2),polVol(:,3));
% er.Color =[0 0 0];
% er.LineStyle = 'none';
% title('Polymer Volensity')
% axis square

figure
hold on
Plotting.shadedErrorBar(xTime,polVol(:,1),polVol(:,2))
title('Polymer Volume')
axis square
box on
xlim([xTime(1) xTime(end)]);
%% Decay comparison
fieldN = fieldnames(data);
idx = 4;
figure
hold on
for i = 1:numel(fieldN)
    currF = ['T' num2str(i-1)];
    plot(data(idx).(currF).intRes{1}.Distance,data(idx).(currF).intRes{1}.normInt);
   
end
legend(xAxis)
box on
xlabel('Distance from cell')
ylabel('Norm. Intensity.')
axis square

%% 

fieldN = fieldnames(data);
figure
hold on
bins = linspace(0,100000,100);
for i = 1:numel(fieldN)
    currF = ['T' num2str(i-1)];
    currAvg = zeros(1,length(bins)-1);
    for j = 1:size(data,2)
       currDist = data(j).(currF).distances;
       [N,edges] = histcounts(currDist{1},bins);
       
       currAvg = currAvg+N;
       
    end
    
    currAvg = currAvg/size(data,2);
    
    plot(edges(2:end)-mean(diff(edges)),N/sum(N));
    
   
    
   
end
legend(xAxis)
box on
xlabel('Distance from cell')
ylabel('Norm. Occurrence')
axis square
ylim([0 0.1])


%% Distance Comparison
fieldN = fieldnames(data);
medDist = zeros(numel(fieldN),3);
maxDist  = zeros(numel(fieldN),3);

for i = 1:numel(fieldN)
    currF = ['T' num2str(i-1)];
    
    tmpMeanDist = zeros(size(data,2),1);
    tmpMaxDist  = zeros(size(data,2),1);
    for j = 1:size(data,2)
       
        tmpMeanDist(j) = median(data(j).(currF).distances{1});
        tmpMaxDist(j)  = max(data(j).(currF).distances{1});
        
    end
    
    medDist(i,1) = mean(tmpMeanDist);
    medDist(i,2) = std(tmpMeanDist);
    medDist(i,3) = std(tmpMeanDist);
    
    maxDist(i,1) = mean(tmpMaxDist);
    maxDist(i,2) = std(tmpMaxDist);
    maxDist(i,3) = std(tmpMaxDist);
    
    
end

figure
subplot(1,2,1)
bar(xAxis,medDist(:,1))
hold on
er = errorbar(xAxis,medDist(:,1),medDist(:,2),medDist(:,3));
er.Color =[0 0 0];
er.LineStyle = 'none';
title('Median Distance')
axis square
ylabel('Distance (nm)')

subplot(1,2,2)
bar(xAxis,maxDist(:,1))
hold on
er = errorbar(xAxis,maxDist(:,1),maxDist(:,2),maxDist(:,3));
er.Color =[0 0 0];
er.LineStyle = 'none';
title('Max distance')
ylabel('Distance (nm)')
axis square




