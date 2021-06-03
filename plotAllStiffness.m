clear;
clc;
close all;

folder2Data = 'D:\Boris_Hongbo\3_different gel conditions';
conditions = {'1K' '1mg' '1K25nm' '5K'};
conditionsVar = {'A','B','C'};

xAxis = categorical(conditions);
%% Load

files = dir(folder2Data);
counter = 1;
data = struct();

c1 = 1;
c2 = 1;
c3 = 1;
c4 = 1;
for i = 1:size(files,1)

    if contains(files(i).name,'Cell') && files(i).isdir
        tmp = load([files(i).folder filesep files(i).name filesep 'results.mat']);
        if contains(files(i).name(end-1:end),conditions{1})
            data(c1).oneK = tmp.res;
            c1 = c1+1;     
        end
        
        if contains(files(i).name,'1K 25nm')
            data(c2).oneK25 = tmp.res;
            c2 = c2+1;            
        end
        
        if contains(files(i).name,'5K')
            data(c3).fiveK = tmp.res;
            c3 = c3+1;            
        end
        
        if contains(files(i).name,'1mg')
            data(c4).oneK1mg = tmp.res;
            c4 = c4+1;            
        end
    end
end


%% Intensity comparison
fieldN = fieldnames(data);
cellInt = zeros(numel(fieldN),3);
polInt  = zeros(numel(fieldN),3);

for i = 1:numel(fieldN)
    currF = fieldN{i};
    
    tmpCellInt = zeros(size(data,2),1);
    tmpPolInt  = zeros(size(data,2),1);
    for j = 1:size(data,2)
       
        tmpCellInt(j) = data(j).(currF).stats.cellInt;
        tmpPolInt(j)  = data(j).(currF).stats.polInt;
        
    end
    
    cellInt(i,1) = mean(tmpCellInt);
    cellInt(i,2) = cellInt(i,1) - min(tmpCellInt);
    cellInt(i,3) = max(tmpCellInt) - cellInt(i,1);
    
    polInt(i,1) = mean(tmpPolInt);
    polInt(i,2) = polInt(i,1) - min(tmpPolInt);
    polInt(i,3) = max(tmpPolInt) - polInt(i,1);
    
   
end

figure
subplot(1,2,1)
bar(xAxis,cellInt(:,1))
hold on
er = errorbar(xAxis,cellInt(:,1),cellInt(:,2),cellInt(:,3));
er.Color =[0 0 0];
er.LineStyle = 'none';
title('Cell intensity')
axis square

subplot(1,2,2)
bar(xAxis,polInt(:,1))
hold on
er = errorbar(xAxis,polInt(:,1),polInt(:,2),polInt(:,3));
er.Color =[0 0 0];
er.LineStyle = 'none';
title('Polymer intensity')
axis square




%% Volume comparison
fieldN = fieldnames(data);
cellVol = zeros(numel(fieldN),3);
polVol  = zeros(numel(fieldN),3);

for i = 1:numel(fieldN)
    currF = fieldN{i};
    
    tmpCellVol = zeros(size(data,2),1);
    tmpPolVol  = zeros(size(data,2),1);
    for j = 1:size(data,2)
       
        tmpCellVol(j) = data(j).(currF).stats.cellVol;
        tmpPolVol(j)  = data(j).(currF).stats.polVol;
        
    end
    
    cellVol(i,1) = mean(tmpCellVol);
    cellVol(i,2) = cellVol(i,1) - min(tmpCellVol);
    cellVol(i,3) = max(tmpCellVol) - cellVol(i,1);
    
    polVol(i,1) = mean(tmpPolVol);
    polVol(i,2) = polVol(i,1) - min(tmpPolVol);
    polVol(i,3) = max(tmpPolVol) - polVol(i,1);
    
   
end

figure
subplot(1,2,1)
bar(xAxis,cellVol(:,1))
hold on
er = errorbar(xAxis,cellVol(:,1),cellVol(:,2),cellVol(:,3));
er.Color =[0 0 0];
er.LineStyle = 'none';
title('Cell Volume')
axis square

subplot(1,2,2)
bar(xAxis,polVol(:,1))
hold on
er = errorbar(xAxis,polVol(:,1),polVol(:,2),polVol(:,3));
er.Color =[0 0 0];
er.LineStyle = 'none';
title('Polymer Volume')
axis square



%% Decay comparison
fieldN = fieldnames(data);
idx = 3;
figure
hold on
for i = 1:numel(fieldN)
    currF = fieldN{i};
    
    plot(data(idx).(currF).intRes{1}.Distance,data(idx).(currF).intRes{1}.normInt);
   
    
   
end
legend(xAxis)
box on
xlabel('Distance from cell')
ylabel('Norm. Intensity.')
axis square
%% Distance Plot

fieldN = fieldnames(data);
idx = 3;
figure
hold on
for i = 1:numel(fieldN)
    currF = fieldN{i};
    
    currDist = data(idx).(currF).distances;
    
    minDist = min(cellfun(@min,currDist));
    maxDist = max(cellfun(@max,currDist));

    bins = linspace(minDist,maxDist,101);
    [N,edges] = histcounts(currDist{1},bins);

    plot(edges(2:end)-mean(diff(edges)),N/sum(N));
    
   
    
   
end
legend(xAxis)
box on
xlabel('Distance from cell')
ylabel('Norm. Occurrence')
axis square
ylim([0 0.05])


%% Distance Comparison
fieldN = fieldnames(data);
medDist = zeros(numel(fieldN),3);
maxDist  = zeros(numel(fieldN),3);

for i = 1:numel(fieldN)
    currF = fieldN{i};
    
    tmpMeanDist = zeros(size(data,2),1);
    tmpMaxDist  = zeros(size(data,2),1);
    for j = 1:size(data,2)
       
        tmpMeanDist(j) = median(data(j).(currF).distances{1});
        tmpMaxDist(j)  = max(data(j).(currF).distances{1});
        
    end
    
    medDist(i,1) = mean(tmpMeanDist);
    medDist(i,2) = medDist(i,1) - min(tmpMeanDist);
    medDist(i,3) = max(tmpMeanDist) - medDist(i,1);
    
    maxDist(i,1) = mean(tmpMaxDist);
    maxDist(i,2) = maxDist(i,1) - min(tmpMaxDist);
    maxDist(i,3) = max(tmpMaxDist) - maxDist(i,1);
    
    
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




