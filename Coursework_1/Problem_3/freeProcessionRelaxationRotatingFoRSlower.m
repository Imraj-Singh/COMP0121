clc
clear
%% Free procession with relaxation in rotating frame of reference, isocromat faster than larmor
% Author: Imraj Singh 03/11/2020

% Given parameters

% Magnitisation aligned along y-axis initially taken from problem 2.4
M = [0, 1, 0]';

% Gyromagnetic ratio for proton in hydrogen
gamma = 2.68*10^8;

% Prescribe a static field strength B0 (T)
B0 = 3;

% Prescribe T1

T2 = 5;

T1 = T2 * 10;

% Calculate neccessary parameters

% Precession frequency
omega0 = 10;%B0*gamma;

% Define modelling parameters

t = 100;%*10^(-3);

% Duration and timestep of model

time = linspace(0,t,1000001);

Msoln = zeros(1000001,3);

Msoln(:,1) = exp(-time./T2).*(M(1)*cos(-omega0/100.*time) + M(2)*sin(-omega0/100.*time));
Msoln(:,2) = exp(-time./T2).*(M(2)*cos(-omega0/100.*time) - M(1)*sin(-omega0/100.*time));
Msoln(:,3) = M(3)*exp(-time./T1) + 1*(1-exp(-time./T1));

plot3(Msoln(:,1),Msoln(:,2),Msoln(:,3))




