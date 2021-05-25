clear;
clc;
close all;

folder2Data = 'D:\Boris_Hongbo\2_after adding cytoD_5 time points\Cell2';
times = [0 2 7 15 24];




%% Loading data

files = dir(folder2Data);
data = struct();
counter = 1;
for i = 1:length(files)
    if ~files(i).isdir && contains(files(i).name,'.mat')
        tmp = load([files(i).folder filesep files(i).name]);
        data.results{counter} = tmp.res;
        counter = counter +1;
    end
    
end

%% Intensity vs Culture time

figure
    subplot(1,2,1)
    plot(times,data.results{1}.stats.cellInt)
    title('Cell intensity')
    xlabel('Culture time')
    ylabel('Avg. Int. (a.u)');
    ylim([0 1.2*max(data.results{1}.stats.cellInt)])
    axis square

    subplot(1,2,2)
    plot(times,data.results{1}.stats.polInt)
    title('Pol. intensity')
    xlabel('Culture time')
    ylabel('Avg. Int. (a.u)');
    ylim([0 1.2*max(data.results{1}.stats.polInt)])
    axis square
 %% Volume vs culture time

%plot
figure
    subplot(1,2,1)
    plot(times,data.results{1}.stats.cellVol)
    title('Cell Volume')
    xlabel('Culture time')
    ylabel('Avg. Int. (a.u)');
    ylim([0 1.2*max(data.results{1}.stats.cellVol)])
    axis square
    
    subplot(1,2,2)
    plot(times,data.results{1}.stats.polVol)
    title('Pol. Volume')
    ylim([0 1.2*max(data.results{1}.stats.polVol)])
    xlabel('Culture time')
    ylabel('Avg. Int. (a.u)');
    axis square
 
 %% Intensity decay vs culture time
 figure
 hold on
 lab = cell(length(times),1);
 for i = 1:length(times) 
    data2Plot = data.results{1}.intRes{i};
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
    data2Plot = data.results{1}.distances(i);
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
 