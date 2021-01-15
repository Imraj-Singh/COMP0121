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
dk = 1/FOV;

% Define the spatial domain over which the data is
x = linspace(-FOV/2, FOV/2 - FOV/n,n);

% Correctly orientated brain
corrected = brain';

% k-space values
k = -n/2*dk:dk:(n/2*dk-dk);

% Define a colour map
cmap = gray(256);
% Do the signal
signal = abs(ifftshift(fft2(fftshift(corrected))));
signal(1,1) = 0;
% Do the image thing
imagesc(k/1000,k/1000,signal)
% Flip the y-axis
set(gca,'YDir','normal')
% Apply the custom colormap.
colormap(cmap); 
% Add colour bar
colorbar;
% Set a log scale
set(gca,'ColorScale','log')
% Make the image aspect ratio square
pbaspect([1 1 1])
% Set the axes
xlabel("$k_x$ (1/mm)", "interpreter", "latex", "fontsize", 10)
ylabel("$k_y$ (1/mm)", "interpreter", "latex", "fontsize", 10)
set(gca, 'defaultAxesTickLabelInterpreter','latex');  