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

% Subplot graph
subplot(2,2,1)
% Image
imagesc(x*1000,x*1000,corrected)
% Flip the y-axis
set(gca,'YDir','normal')
% Make the image aspect ratio square
pbaspect([1 1 1])
% Set the axes
xlabel("$x$ (mm)", "interpreter", "latex", "fontsize", 10)
ylabel("$y$ (mm)", "interpreter", "latex", "fontsize", 10)
set(gca, 'defaultAxesTickLabelInterpreter','latex'); 

% Get k-space values
kspace = ifftshift(fft2(fftshift(corrected)));
plotting = kspace;
plotting(1,1) = 0;
% Subplot graph
subplot(2,2,2)
% Image
imagesc(k/1000,k/1000,abs(plotting))
% Flip the y-axis
set(gca,'YDir','normal')
% Apply the custom colormap.
colormap(cmap); 
% Set a log scale
set(gca,'ColorScale','log')
% Make the image aspect ratio square
pbaspect([1 1 1])
% Set the axes
xlabel("$k_x$ (1/mm)", "interpreter", "latex", "fontsize", 10)
ylabel("$k_y$ (1/mm)", "interpreter", "latex", "fontsize", 10)
set(gca, 'defaultAxesTickLabelInterpreter','latex');  

% Create the rect filter
rect = zeros(n,n);
W = 32*dk;
rect(k > W/2,:)=1;
rect(k < -W/2,:)=1;
rect(:,k > W/2)=1;
rect(:,k < -W/2)=1;

% Apply the rect function
highpasskspace = rect.*kspace;
plotting = highpasskspace;
plotting(1,1) = 0;
plotting(1,2) = max(abs(kspace(:)));
% Subplot graph
subplot(2,2,3)
% Image
imagesc(k/1000,k/1000,abs(plotting))
% Flip the y-axis
set(gca,'YDir','normal')
% Apply the custom colormap.
colormap(cmap); 
% Set a log scale
set(gca,'ColorScale','log')
% Make the image aspect ratio square
pbaspect([1 1 1])
% Set the axes
xlabel("$k_x$ (1/mm)", "interpreter", "latex", "fontsize", 10)
ylabel("$k_y$ (1/mm)", "interpreter", "latex", "fontsize", 10)
set(gca, 'defaultAxesTickLabelInterpreter','latex');  

% Create the image of the effect of the filter
highpassimage = fftshift(ifft2(ifftshift(highpasskspace)));
% Subplot graph
subplot(2,2,4)
% Image
imagesc(x*1000,x*1000,abs(highpassimage))
% Flip the y-axis
set(gca,'YDir','normal')
% Apply the custom colormap.
colormap(cmap); 
% Make the image aspect ratio square
pbaspect([1 1 1])
% Set the axes
xlabel("$x$ (mm)", "interpreter", "latex", "fontsize", 10)
ylabel("$y$ (mm)", "interpreter", "latex", "fontsize", 10)
set(gca, 'defaultAxesTickLabelInterpreter','latex'); 

