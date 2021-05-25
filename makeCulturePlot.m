clear;
clc;
close all;

folder2Data = 'D:\Boris_Hongbo\1_with culturing time_0-24h\Cell 2 - Results';
times = [0 2 7 15 24];




%% Loading data

files = dir(folder2Data);
data = struct();
counter = 1;
for i = 1:length(files)
    if ~files(i).isdir && ~contains(files(i).name,'.fig')
        tmp = load([files(i).folder filesep files(i).name]);
        data.results{counter} = tmp.res;
        counter = counter +1;
    end
    
end

assert(length(times) == length(data.results),'Number of time points does not match number of files in directory');

%% Intensity vs Culture time
%extract data
cellData = zeros(1,length(times));
polData  = zeros(1,length(times));
for i = 1:length(times)
    cellData(i) = data.results{i}.stats.cellInt; 
    polData(i)  = data.results{i}.stats.polInt;
    
end

figure
    subplot(1,2,1)
    plot(times,cellData)
    title('Cell intensity')
    xlabel('Culture time')
    ylabel('Avg. Int. (a.u)');
    axis square

    subplot(1,2,2)
    plot(times,polData)
    title('Pol. intensity')
    xlabel('Culture time')
    ylabel('Avg. Int. (a.u)');
    axis square
 %% Volume vs culture time
%extract data
cellData = zeros(1,length(times));
polData  = zeros(1,length(times));
for i = 1:length(times)
    cellData(i) = data.results{i}.stats.cellVol; 
    polData(i)  = data.results{i}.stats.polVol;
    
end
%plot
figure
    subplot(1,2,1)
    plot(times,cellData)
    title('Cell Volume')
    xlabel('Culture time')
    ylabel('Avg. Int. (a.u)');
    axis square
    
    subplot(1,2,2)
    plot(times,polData)
    title('Pol. Volume')
    xlabel('Culture time')
    ylabel('Avg. Int. (a.u)');
    axis square
 
 %% Intensity decay vs culture time
 figure
 hold on
 lab = cell(length(times),1);
 for i = 1:length(times) 
    data2Plot = data.results{i}.intRes;
    plot(data2Plot.Distance,data2Plot.normInt)
    lab{i} = [num2str(times(i)) 'h']; 
end
 title('Intensity decay')
 xlabel('Distance from the cell (\mum)')
 ylabel('Norm. Intensity (a.u.)')
 legend(lab);
 box on;
 axis square
 
 %% Distance distribution vs culture time
 
 figure
 hold on
 lab = cell(length(times),1);
 for i = 1:length(times) 
    data2Plot = data.results{i}.distances;
    minDist = min(cellfun(@min,data2Plot));
    maxDist = max(cellfun(@max,data2Plot));

    bins = linspace(minDist,maxDist,101);
    [N,edges] = histcounts(data2Plot{1},bins);
    
    plot(edges(2:end)-mean(diff(edges)),N/sum(N));
    
    lab{i} = [num2str(times(i)) 'h']; 
end
 title('Densification range')
 xlabel('Distance (\mum)')
 ylabel('Normalized distribution')
 legend(lab);
 box on;
 axis square
 