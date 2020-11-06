clc
clear

%% animation of spin precession
% Initialize the spin vector
% magnitude in the unit of gamma*h_bar
mu = 1/2;
% set the initial direction in terms of polar and azimuthal angles
% polar angle
theta = pi/6;
% azimuthal angle
phi = 0;
% compute the Cartesian components of the vector
vecMu = mu*[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)]';

% set the precession frequency relative to the Larmor frequency
% in the unit of radian*kHz
omega = 2*pi;
% set the time increment (in the unit of ms)
deltaT = 0.01;
% set the number of time increments to animate
noOfSteps = 200;
% amount rotated in a time increment
deltaphi = 0.1;

%% Ploting

h1 = figure;
hold on;
axis equal;
view(100, 10);
xlabel('\mu_x');
ylabel('\mu_y');
zlabel('\mu_z');
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
zlim([0 0.5]);
grid on;

for i=0:noOfSteps-1
    plot3([0 vecMu(1)], [0 vecMu(2)], [0 vecMu(3)], 'c-','LineWidth',2);
    plot3([0 vecMu(1)], [0 vecMu(2)], [0 vecMu(3)], 'r.','MarkerSize',10);
    pause(0.05)
    vecMu = rotateZ(vecMu,deltaphi);
end

