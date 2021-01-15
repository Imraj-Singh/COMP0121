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

% Define Gaussian
gaus = @(x,mu,sig)(1/(sig*sqrt(2)))*exp(-(((x-mu).^2)/(2*sig.^2)));

% Define hann
hann = [zeros(1,nnz(k<-8*dk)),(0.5 - 0.5*cos(2*pi*(1:17)/18)),zeros(1,nnz(k>8*dk))];

% Calculate signal
signal = gaus(k,0,8*dk).*hann;


% original
h = subplot(3,1,1);
plotComplex(k/1000,signal,h);
grid on
box on
hold off
title('Signal = Gaussian * Hanning', "interpreter", "latex", "fontsize", 10)
xlabel("$k$ (1/mm)", "interpreter", "latex", "fontsize", 10)
ylabel("Signal", "interpreter", "latex", "fontsize", 10)

% brute force
h = subplot(3,1,2);
plotComplex(x*1000,fftshift(ifft(ifftshift(hann))*N),h);
grid on
box on
hold off
title("Point spread function of Hanning function", "interpreter", "latex", "fontsize", 10)
% Set the axes labels
xlabel("$x$ (mm)", "interpreter", "latex", "fontsize", 10)
ylabel("$\rho$", "interpreter", "latex", "fontsize", 10)

h = subplot(3,1,3);
plotComplex(x*1000,fftshift(ifft(ifftshift(signal))*N),h);
grid on
box on
hold off
title("IFFT (Signal)", "interpreter", "latex", "fontsize", 10)
% Set the axes labels
xlabel("$x$ (mm)", "interpreter", "latex", "fontsize", 10)
ylabel("$\rho$", "interpreter", "latex", "fontsize", 10)

