clc
clear
%% Forced procession ignoring relaxation
% Author: Imraj Singh 03/11/2020

% Given parameters

% Magnitisation aligned along z-axis initially
M = [0, 0, 1]';

% Duration of the RF pulse
t = 1%*10^(-3);

% Flip angle
dtheta = 90*pi/180;

% Gyromagnetic ratio for proton in hydrogen
gamma = 2.68*10^8;

% Calculate neccessary parameters

% The magnetic field strength caused by RF
B1 = dtheta/(gamma*t);

% Precession frequency
omega1 = gamma*B1;

% Define modelling parameters

% Duration and timestep of model

time = linspace(0,t,101);

Msoln = zeros(101,3);

Msoln(:,1) = M(1);
Msoln(:,2) = M(2)*cos(omega1.*time) + M(3)*sin(omega1.*time);
Msoln(:,3) = M(3)*cos(omega1.*time) - M(2)*sin(omega1.*time);

plot3(Msoln(:,1),Msoln(:,2),Msoln(:,3))




