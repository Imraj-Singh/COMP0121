clc
clear
%% Free procession ignoring relaxation
% Author: Imraj Singh 03/11/2020

% Given parameters

% Magnitisation aligned along y-axis initially taken from problem 2.4
M = [0, 1, 0]';

% Gyromagnetic ratio for proton in hydrogen
gamma = 2.68*10^8;

% Prescribe a static field strength B0 (T)
B0 = 3;

% Calculate neccessary parameters

% Precession frequency
omega0 = B0*gamma;

% Define modelling parameters

t = 1*10^(-3);

% Duration and timestep of model

time = linspace(0,t,1001);

Msoln = zeros(1001,3);

Msoln(:,1) = M(1).*(M(1)*cos(omega0.*time) + M(2)*sin(omega0.*time));
Msoln(:,2) = M(2)*cos(omega0.*time) - M(1)*sin(omega0.*time);
Msoln(:,3) = 0;%M(3)*cos(omega1.*time) - M(2)*sin(omega1.*time);

plot3(Msoln(:,1),Msoln(:,2),Msoln(:,3))




