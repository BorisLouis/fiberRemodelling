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
    cellInt(i,:) = data.(currF).stats.cellInt(1:nCond)-data.(currF).stats.cellCornerInt(1:nCond);
    polInt(i,:)  = data.(currF).stats.polInt(1:nCond)./data.(currF).stats.polCornerInt(5);
    
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
    cellVol(i,:) = cellVol(i,:)./cellVol(i,1);
    polVol(i,:)  = data.(currF).stats.polVol(1:nCond);
    polVol(i,:)  = polVol(i,:)./polVol(i,1);
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

%% Volume normalized
fieldN = fieldnames(data);
nData = numel(fieldN);
nCond = length(xTime);
cellVol = zeros(nData,nCond);
polVol  = zeros(nData,nCond);
normPolVol = polVol;
for i = 1:numel(fieldN)
    currF = ['C' num2str(i)];
    cellVol(i,:) = data.(currF).stats.cellVol(1:nCond);
    polVol(i,:)  = data.(currF).stats.polVol(1:nCond);
    normPolVol(i,:) = polVol(i,:)./cellVol(i,:);
    
end

% shadow error bar
figure
Plotting.shadedErrorBar(xTime,median(normPolVol,1),std(normPolVol,1))
title('Normalized polymer volume')
xlim([xTime(1) xTime(end)]);
xlabel('Time (min)')
ylabel('Pol. Vol. Rel. to cell')
axis square
box on
%% Decay comparison
fieldN = fieldnames(data);
idx = 1;
figure
hold on
for i = 1:numel(fieldN)
    currF = ['C' num2str(i)];
    plot(data.(currF).intRes{idx}.Distance,data.(currF).intRes{idx}.normInt);
   
end
legend(xAxis)
box on
xlabel('Distance from cell')
ylabel('Norm. Intensity.')
axis square



% bin distance curve to average it out
cellDecay = cell(2,numel(fieldN));
for i = 1:size(data.C1.intRes,1)
    for j = 1:numel(fieldN)
        currF = ['C' num2str(j)];
        cellDecay{1,i} = [cellDecay{1,i}; data.(currF).intRes{i}.Distance(:)];
        cellDecay{2,i} = [cellDecay{2,i}; data.(currF).intRes{i}.normInt(:)];
        
        
    end
    
end

bins = linspace(0,400000,1000);

for i = 1:size(data.C1.intRes,1)
   % sort x axis
   [cellDecay{1,i}, Ia] = sort(cellDecay{1,i});
   cellDecay{2,i} = cellDecay{2,i}(Ia); 
    
    [~,~,loc]=histcounts(cellDecay{1,i},bins);
    meanInt  = accumarray(loc(loc>0),cellDecay{2,i})./accumarray(loc(loc>0),1);
    cellDecay{3,i} = unique(bins(loc(loc>0)));
    cellDecay{4,i} = meanInt(~isnan(meanInt));    
end




figure
hold on
for i = 1:numel(fieldN)
   plot(cellDecay{3,i},cellDecay{4,i})
   
end
legend(xAxis)
box on
xlabel('Distance from cell')
ylabel('Norm. Intensity.')
axis square


%% Distance distribution
fieldN = fieldnames(data);
figure
hold on
bins = linspace(0,100000,100);

for i = 1:size(data.C1.intRes,1)
    currAvg = zeros(1,length(bins)-1);
    for j = 1:numel(fieldN)
       currF = ['C' num2str(j)];
       currDist = data.(currF).distances{i};
       [N,edges] = histcounts(nonzeros(currDist),bins);
       
       currAvg = currAvg+N;
       
    end
    
    currAvg = currAvg/sum(currAvg);
    
    plot(edges(2:end)-(mean(diff(edges)/2)),currAvg);
    

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
for i = 1:size(data.C1.intRes,1)
    
    tmpMeanDist = zeros(numel(fieldN),1);
    tmpMaxDist  = zeros(numel(fieldN),1);
    for j = 1:numel(fieldN)
        currF = ['C' num2str(j)];
        tmpMeanDist(j) = median(nonzeros(data.(currF).distances{i}));
        tmpMaxDist(j)  = max(data.(currF).distances{i});
        
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

%% Distance T0-T5
bins = linspace(0,100000,100);
figure
hold on
for i = 1:numel(fieldN)
   currF = ['C' num2str(i)];
   currAvg = zeros(1,length(bins)-1);
   currDist = data.(currF).polDist;
   [N,edges] = histcounts(currDist,bins);

   plot(edges(1:end-1),N)
end
legend({'C1','C2','C3','C4','C5'})
xlabel('Distance (nm)')
ylabel('Occurences')
set(gca,'XScale','log')
set(gca,'YScale','log');
axis square
box on


%% Model T0-T5
idx = 4;
currF = ['C' num2str(idx)];
data2Render = data.(currF).cellMask(:,:,:,1);
iSurface = isosurface(data2Render,1/2);

figure(idx)
hold on
p2 = patch(iSurface);
p2.FaceColor = [0,1,0];
p2.EdgeColor = 'none';

clear data2Render iSurface;

data2Render = data.(currF).polymerMask(:,:,:,1);
iSurface = isosurface(data2Render,1/2);

p3 = patch(iSurface);
p3.FaceColor = [1 0 0];
p3.FaceAlpha = 0.2;
p3.EdgeColor = 'none';

clear data2Render iSurface;

data2Render = data.(currF).polymerMask(:,:,:,end);
iSurface = isosurface(data2Render,1/2);

p3 = patch(iSurface);
p3.FaceColor = [0.7 0.7 0.7];
p3.FaceAlpha = 0.5;
p3.EdgeColor = 'none';

view(3);
axis tight
camlight
lighting gouraud
title('unicolor');




