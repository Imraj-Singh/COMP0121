clc
clear
close all
% author: Imraj Singh 06/02/2021

% Add path to data
addpath 'C:\Users\imraj\Documents\UCL\COMP0121\Coursework_4\Brain_Reconstruction'

% load data
load('data.mat')

n = 54
N_k = 2000^2;                         % number of k-samples
t   = linspace(0,sqrt(0.5),N_k)';   % dummy variable to parameterise spiral
k_x = (1000).*t.^2.*cos(2*pi*n*t);  % spiral kx-coords
k_y = (1000).*t.^2.*sin(2*pi*n*t);  % spiral ky-coords

% set the parametes
dx = 1/1000;
n = 256;
FOV = dx*n;
dk = 1/FOV;

k_x = floor(k_x/dk);
k_y = floor(k_y/dk);
k_x = k_x + 129;
k_y = k_y + 129;


%%

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
subplot(1,2,1)
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
% Subplot graph
subplot(1,2,2)
% Image
imagesc(k/1000,k/1000,abs(kspace))
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



%%

% Reject points under a certain values

check = zeros('like',kspace);
for i=1:N_k
     check(k_x(i),k_y(i)) = kspace(k_x(i),k_y(i));
end
kspace = check;


% Subplot graph
subplot(1,2,1)
% Image
imagesc(k/1000,k/1000,abs(kspace))
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
highpassimage = fftshift(ifft2(ifftshift(kspace)));
% Subplot graph
subplot(1,2,2)
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


% Requires R2020a or later
exportgraphics(gcf,'myplot.png','Resolution',300) 

1-length(nonzeros(kspace))/(256^2)
