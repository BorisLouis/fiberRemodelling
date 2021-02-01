clear;
close all;
clc;
%% User Input
file.path = 'D:\Documents\Unif\PhD\2020-Data\07- July\MAX\Example data for Boris_1k_10nm_TFM\T4';
file.ext  = '.tif';

info.pxSizeXY = 568; 
info.pxSizeZ  = 569;

%Give info about the channels, the word needs to be lowercase with no typos
chan.ch00 = 'cell';
chan.ch01 = 'polymer';
chan.ch02 = 'nucleus';
chan.ch03 = 'ignore';

%% Loading data
stack = Core.fiberRemodelling(file,info);

stack.loadData(chan);

stack.showChannel;

%% analysis occurs here

Mask = stack.calc3DMask();
%%
stack.plotCellContour();

%% 3D
%this step takes time so we don't really run it
%stack.renderCell3D();

%% Intensity analysis
weight = [info.pxSizeXY,info.pxSizeXY,info.pxSizeZ];

step = 2*info.pxSizeXY;
%test = DistMap.calcWeightedDistMap(Mask,weight);
stack.intensityDistrib(weight,step);

