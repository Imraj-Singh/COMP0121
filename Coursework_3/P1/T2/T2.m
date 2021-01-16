clc
clear
close all
% author: Imraj Singh 06/01/2021

% Define the Field of view for the simulation
FOV = 8/1000;

% Calculate the k-space spacing
dk = 1/FOV;

% Define the number of points according to 2^p rule needed for IFFT
p = 7;
N = 2^p;

% Define the frequencies
k = linspace(-dk*(N/2),dk*(N/2-1),N);

% Define the distances
x = linspace(-FOV/2,FOV/2 - FOV/N,N);

% Impulse location
val = dk*0;

% Make signal vector with unit impulse
signal = double(k == val);

figure
% original
h = subplot(3,1,1);
plotComplex(k/1000,signal,h);
grid on
box on
hold off
title('Signal', "interpreter", "latex", "fontsize", 10)
xlabel("$k$ (1/mm)", "interpreter", "latex", "fontsize", 10)
ylabel("Signal", "interpreter", "latex", "fontsize", 10)

% original
h = subplot(3,1,2);
plotComplex(x*1000,ifft(signal),h);
grid on
box on
hold off
ylim([-0.01 0.01])
title('IFFT(Signal) - Not shifted', "interpreter", "latex", "fontsize", 10)
xlabel("$x$ (mm)", "interpreter", "latex", "fontsize", 10)
ylabel("$\rho$", "interpreter", "latex", "fontsize", 10)

% shifted
h = subplot(3,1,3);
plotComplex(x*1000,fftshift(ifft(ifftshift(signal))),h);
grid on
box on
hold off
ylim([-0.01 0.01])
title('IFFT(Signal) - Shifted', "interpreter", "latex", "fontsize", 10)
% Set the axes labels
xlabel("$x$ (mm)", "interpreter", "latex", "fontsize", 10)
ylabel("$\rho$", "interpreter", "latex", "fontsize", 10)

