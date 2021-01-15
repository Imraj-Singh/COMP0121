clc
clear
close all
% author: Imraj Singh 06/01/2021

% Add path to data
addpath 'C:\Users\Imraj Singh\Documents\UCL\Comp_MRI\COMP0121\Coursework_3\data'

% load data
load('data.mat')

% set the parametes
dx = 1/1000;
n = 256;
FOV = dx*n;

% Define the spatial domain over which the data is
x = linspace(-FOV/2, FOV/2 - FOV/n,n);

% Define a colour map
cmap = gray(256);
% Flip and plot brain data so that it is in the correct orientation
imagesc(x*1000,x*1000,brain')
% Flip the y-axis
set(gca,'YDir','normal')
% Apply the custom colormap.
colormap(cmap);
% Add colour bar
colorbar;
% Make the image aspect ratio square
pbaspect([1 1 1])
% Set the axes
xlabel("$x$ (mm)", "interpreter", "latex", "fontsize", 10)
ylabel("$y$ (mm)", "interpreter", "latex", "fontsize", 10)
set(gca, 'defaultAxesTickLabelInterpreter','latex');  