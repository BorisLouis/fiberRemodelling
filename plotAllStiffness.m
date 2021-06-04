clear;
clc;
close all;

folder2Data = 'D:\Boris_Hongbo\3_different gel conditions';
conditions = {'1K1mg' '1K1.5mg' '1K25nm' '5K'};

xAxis = categorical(conditions,'Ordinal',1);
xAxis = reordercats(xAxis,conditions);
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
        if contains(files(i).name(end-1:end),'1K')
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

fieldN = fieldnames(data);
fieldN = fieldN([3 1 2 4]);


%% Intensity comparison

cellInt = zeros(numel(fieldN),3);
polInt  = zeros(numel(fieldN),3);

for i = 1:numel(fieldN)
    currF = fieldN{i};
    idx = ~cellfun(@isempty,{data.(currF)});
    nData = sum(idx);
    tmpCellInt = zeros(nData,1);
    tmpPolInt  = zeros(nData,1);
    for j = 1:nData
       
        tmpCellInt(j) = data(j).(currF).stats.cellInt;
        tmpPolInt(j)  = data(j).(currF).stats.polInt;
        
    end
    
    cellInt(i,1) = mean(tmpCellInt);
    cellInt(i,2) = std(tmpCellInt);
    cellInt(i,3) = std(tmpCellInt);
    
    polInt(i,1) = mean(tmpPolInt);
    polInt(i,2) = std(tmpPolInt);
    polInt(i,3) = std(tmpPolInt);
    
   
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
cellVol = zeros(numel(fieldN),3);
polVol  = zeros(numel(fieldN),3);
normPolVol = polVol;
for i = 1:numel(fieldN)
    currF = fieldN{i};
    idx = ~cellfun(@isempty,{data.(currF)});
    nData = sum(idx);
    tmpCellVol = zeros(nData,1);
    tmpPolVol  = zeros(nData,1);
    tmpNormPolVol = zeros(nData,1);
    
    for j = 1:nData
       
        tmpCellVol(j) = data(j).(currF).stats.cellVol;
        tmpPolVol(j)  = data(j).(currF).stats.polVol;
        tmpNormPolVol(j) = data(j).(currF).stats.polVol./data(j).(currF).stats.cellVol;
    end
    
    cellVol(i,1) = mean(tmpCellVol);
    cellVol(i,2) = std(tmpCellVol);
    cellVol(i,3) = std(tmpCellVol);
    
    polVol(i,1) = mean(tmpPolVol);
    polVol(i,2) = std(tmpPolVol);
    polVol(i,3) = std(tmpPolVol);
    
    normPolVol(i,1) = mean(tmpNormPolVol);
    normPolVol(i,2) = std(tmpNormPolVol);
    normPolVol(i,3) = std(tmpNormPolVol);
      
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


%normalized Volume
figure
bar(xAxis,normPolVol(:,1))
hold on
er = errorbar(xAxis,normPolVol(:,1),normPolVol(:,2),normPolVol(:,3));
er.Color =[0 0 0];
er.LineStyle = 'none';
title('Polymer Volume')
axis square

%% Decay comparison
% bin distance curve to average it out

cellDecay = cell(2,numel(fieldN));
for i = 1:numel(fieldN)
    currF = fieldN{i};
    
    idx = ~cellfun(@isempty,{data.(currF)});
    nData = sum(idx);
    
    for j = 1:nData
        
        cellDecay{1,i} = [cellDecay{1,i}; data(j).(currF).intRes{1}.Distance(:)];
        cellDecay{2,i} = [cellDecay{2,i}; data(j).(currF).intRes{1}.normInt(:)];
        
    end
end

bins = linspace(0,400000,1000);

for i = 1:size(cellDecay,2)
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
figure
hold on
bins = linspace(0,100000,100);

for i = 1:numel(fieldN)
    currAvg = zeros(1,length(bins)-1);
    currF = fieldN{i};
    idx = ~cellfun(@isempty,{data.(currF)});
    nData = sum(idx);
    for j = 1:nData
       currDist = data(j).(currF).distances{1};
       [N,edges] = histcounts(nonzeros(currDist),bins);
       
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

medDist = zeros(numel(fieldN),3);
maxDist  = zeros(numel(fieldN),3);

for i = 1:numel(fieldN)
    currF = fieldN{i};
    idx = ~cellfun(@isempty,{data.(currF)});
    nData = sum(idx);
    tmpMeanDist = zeros(nData,1);
    tmpMaxDist  = zeros(nData,1);
    for j = 1:nData
       
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




