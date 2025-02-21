clear;
close all;
clc;
%% User Input
file.path = 'K:\different peptides\test';
%file.path = [path filesep 'Split'];
file.ext  = '';

info.pxSizeXY = 570; 
info.pxSizeZ  = 570;

%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.c001 = 'cell';
chan.c002 = 'polymer';
chan.c003 = '';

mkdir( [file.path filesep 'Figures']);
%% Loading data


stack = Core.fiberRemodelling(file,info);

stack.loadData(chan);

stack.showChannel;

%% analysis occurs here

Mask = stack.calc3DMask();

%%
%stack.plotCellContour();

%% 3D
%this step takes time so we don't really run it
try
    stack.renderCell3D(1);
catch
    rendering3D.compile3DRendering();
    stack.renderCell3D(1);
end 
%% Gel network

stack.getDensifiedNetwork();


%% Display polymer and cell in 3D
stack.renderCellPolymer3D(1)



%% Get cell 'hole' in polymer

stack.getCellGap();

%% display all
stack.renderAll()


%% 

[Stats] = stack.calcStats();


%% Intensity analysis
% weight = [info.pxSizeXY,info.pxSizeXY,info.pxSizeZ];
% 
% step = 2*info.pxSizeXY;
% %test = DistMap.calcWeightedDistMap(Mask,weight);
% stack.intensityDistrib(weight,step);

%% get distance between polymer and cell

distance = stack.calcDistCellPol;


%% Plot Distance
minDist = min(cellfun(@min,distance));
maxDist = max(cellfun(@max,distance));

bins = linspace(minDist,maxDist,101);
leg=cell(length(distance),1);
figure
hold on
for i = 1:length(distance)
    [N,edges] = histcounts(distance{i},bins);
    
    bar(edges(2:end)-mean(diff(edges)),N/sum(N));
    
    leg{i} = ['C' num2str(i)];
    
end
xlabel('Distance (\mum)')
ylabel('Area normalized distribution');
xlim([info.pxSizeXY  max(edges)])
axis square
box on
legend(leg)

%% Save data
fileName = [file.path filesep 'results.mat'];
res = stack.getResults();
save(fileName,'res');
f = msgbox('Data succesfully saved');


