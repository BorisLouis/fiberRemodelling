clear;
close all;
clc;
%% User Input
file.path = 'D:\Documents\Unif\PhD\2021-Data\01 - Jan\28 - Hongbo\Boris_Hongbo\cell2\SplitData';
file.ext  = '.tif';

info.pxSizeXY = 568; 
info.pxSizeZ  = 569;

%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.c001 = 'polymer';
chan.c002 = 'cell';
chan.c003 = 'ignore';
chan.c004 = 'ignore';

%% Loading data
rendering3D.compile3DRendering();

stack = Core.fiberRemodelling(file,info);

stack.loadData(chan);

stack.showChannel;

%% analysis occurs here

Mask = stack.calc3DMask();

%%
stack.plotCellContour();

%% 3D
%this step takes time so we don't really run it
stack.renderCell3D(1);

%% Gel network

stack.getDensifiedNetwork();


%% Display polymer and cell in 3D
stack.renderCellPolymer3D(1)


%% 

[Volumes] = stack.calcVolumes();


%%
[NOP]     = stack.calcNOP();


%% Intensity analysis
weight = [info.pxSizeXY,info.pxSizeXY,info.pxSizeZ];

step = 2*info.pxSizeXY;
%test = DistMap.calcWeightedDistMap(Mask,weight);
stack.intensityDistrib(weight,step);

