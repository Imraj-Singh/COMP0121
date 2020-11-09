clc
clear

%% Problem 4.1 FID sequence signal visualisation
% Author: Imraj Singh 03/11/2020

% Given parameters

% Magnitisation aligned along z-axis initially
M = [0, 0, 1]';

% Gyromagnetic ratio for proton in hydrogen
gamma = 2.68*10^8;

% Prescribe T1
T2 = 1;
T1 = T2 * 2;

% Flip angle
dtheta = 90 * pi / 180;

% Calculate neccessary parameters

% Precession frequency
omega0 = 2 * 2 * pi;

s0 = 1;
h = figure;

t=0:0.01:10;

plot(t,sin(omega0.*t).*s0.*exp(-t/T2),'k-','linewidth',2)
hold on
title('Time from $\pi / 2$ flip: ', "interpreter", "latex", "fontsize", 15)
plot(t,s0.*exp(-t/T2),'k-','linewidth',2)

grid on
box on
xlim([0 max(t)]);
ylim([-1 1]);
xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
ylabel("Signal", "interpreter", "latex", "fontsize", 15)
hold off
