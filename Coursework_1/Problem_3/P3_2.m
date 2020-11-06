clc
clear
%% Free procession with relaxation
% Author: Imraj Singh 03/11/2020

% Given parameters

% Magnitisation aligned along y-axis initially taken from problem 2.4
M = [0, 1, 0]';

% Gyromagnetic ratio for proton in hydrogen
gamma = 2.68*10^8;

% Prescribe a static field strength B0 (T)
B0 = 3;

% Prescribe T1

T2 = .5;

T1 = T2 * 2;

% Calculate neccessary parameters

% Precession frequency
omega0 = 10;%B0*gamma;

% Define modelling parameters

t = 5;%*10^(-3);

% Duration and timestep of model

time = linspace(0,t,101);

Msoln = zeros(length(time),3);

Msoln(:,1) = exp(-time./T2).*(M(1)*cos(omega0.*time) + M(2)*sin(omega0.*time));
Msoln(:,2) = exp(-time./T2).*(M(2)*cos(omega0.*time) - M(1)*sin(omega0.*time));
Msoln(:,3) = M(3)*exp(-time./T1) + 1*(1-exp(-time./T1));

%% Animation module

video = VideoWriter(['3_2', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = 100;
video.set('FrameRate', frameRate);

video.open();

h = figure;
frame = getframe(h);
video.writeVideo(frame);

for i=1:length(time)
    quiver3(0,0,0,Msoln(i,1),Msoln(i,2),Msoln(i,3),'linewidth',5,'LineStyle','-','Color','r')
    title(['Time: ', num2str(time(i)),'s'], "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([0 1]);
    xlabel("$M_x$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_y$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_z$", "interpreter", "latex", "fontsize", 15)
    hold on
    plot3(Msoln(1:i-1,1),Msoln(1:i-1,2),Msoln(1:i-1,3),'b-')
    %quiver3(zeros(i-1,1),zeros(i-1,1),zeros(i-1,1),Msoln(1:i-1,1),Msoln(1:i-1,2),Msoln(1:i-1,3),'linewidth',1,'LineStyle','--','Color','b')
    hold off
    frame = getframe(h);
    video.writeVideo(frame);
end

video.close();




