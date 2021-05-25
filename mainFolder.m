clear;
close all;
clc;

path = 'D:\Boris_Hongbo\2_after adding cytoD_5 time points';
info.pxSizeXY = 570; 
info.pxSizeZ  = 570;

%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.c001 = 'cell';
chan.c002 = 'polymer';
chan.c003 = 'ignore';
chan.c004 = 'ignore';


%% get folder to Analyze
folders = dir(path);
idx =[];
for i=1:size(folders,1)
   currentFolder = [folders(i).folder filesep folders(i).name];
   
   testPath = [currentFolder filesep 'Split'];
   
   if isfolder(testPath)
       idx = [idx i];
       
   end
   
end

folder2Analyze = folders(idx);
rendering3D.compile3DRendering();
%% data processing
for i = 1:size(folder2Analyze,1)
    currentFolder = [folder2Analyze(i).folder filesep folder2Analyze(i).name];
    file.path = [currentFolder filesep 'Split'];
    file.ext  = '.tif';

    %% Loading data
    

    stack = Core.fiberRemodelling(file,info);

    stack.loadData(chan);

    %% analysis occurs here

    Mask = stack.calc3DMask();

    %% Gel network

    stack.getDensifiedNetwork();

    %% rendering
    stack.renderCellPolymer3D(1)
    fileName = [currentFolder filesep 'rendering.png'];
    saveas(gcf,fileName)
    
    fileName = [currentFolder filesep 'rendering.fig'];
    saveas(gcf,fileName)
    
    
    
    %% 
    [Stats] = stack.calcStats();

    %% Intensity analysis
    weight = [info.pxSizeXY,info.pxSizeXY,info.pxSizeZ];

    step = 2*info.pxSizeXY;
    %test = DistMap.calcWeightedDistMap(Mask,weight);
    stack.intensityDistrib(weight,step);

    %% get distance between polymer and cell

    distance = stack.calcDistCellPol;

    %% Save data
    fileName = [currentFolder filesep 'results.mat'];
    res = stack.getResults();
    save(fileName,'res');
    close all;
end
