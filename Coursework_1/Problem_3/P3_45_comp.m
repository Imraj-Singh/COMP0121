% housekeeping
clc
clear
%% Free procession with relaxation rotating FoR at Larmor
% Faster than Larmor by omega0/10
% Author: Imraj Singh 03/11/2020

% Given parameters
% Magnitisation aligned along y-axis initially
M = [0, 1, 0]';
M0 = [0, 0, 1]';

% Prescribe T1
T2 = 1;
T1 = T2 * 2;

% Calculate neccessary parameters

% Precession frequency
omega0 = 2 * 2 * pi;

% End time
t = 6;
resolution = 50;

% Duration and timestep of model
time = linspace(0, t, t * resolution + 1);

Msolnf(:,1) = exp(-time./T2).*(M(1)*cos(omega0*11/10.*time) + M(2)*sin(omega0*11/10.*time));
Msolnf(:,2) = exp(-time./T2).*(M(2)*cos(omega0*11/10.*time) - M(1)*sin(omega0*11/10.*time));
Msolnf(:,3) = M(3)*exp(-time./T1) + M0(3)*(1-exp(-time./T1));

Msolns(:,1) = exp(-time./T2).*(M(1)*cos(omega0*9/10.*time) + M(2)*sin(omega0*9/10.*time));
Msolns(:,2) = exp(-time./T2).*(M(2)*cos(omega0*9/10.*time) - M(1)*sin(omega0*9/10.*time));
Msolns(:,3) = M(3)*exp(-time./T1) + M0(3)*(1-exp(-time./T1));




%% Animation module - just the 3D
% Author: Imraj Singh 03/11/2020

% initialise the video
video = VideoWriter(['3_45_comp', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = resolution;
video.set('FrameRate', frameRate);

video.open();

% initialise the figure
h = figure;

% specify animation captures each degree of rotation
for i=1:length(time)
    quiver3(0,0,0,Msolns(i,1),Msolns(i,2),Msolns(i,3),'linewidth',4,'LineStyle','-','Color','r')
    hold on
    quiver3(0,0,0,Msolnf(i,1),Msolnf(i,2),Msolnf(i,3),'linewidth',4,'LineStyle','-','Color','b')
    quiver3(0,0,0,0,1,0,'linewidth',2,'LineStyle','--','Color','k')
    plot3(Msolns(1:i,1),Msolns(1:i,2),Msolns(1:i,3),'r--','linewidth',2)
    plot3(Msolnf(1:i,1),Msolnf(1:i,2),Msolnf(1:i,3),'b--','linewidth',2)

    % format title, legend, labels, limits, grid and box
    title(['Time: ', num2str(time(i)),' s'], "interpreter", "latex", "fontsize", 15)
    legend('Slower','Faster', "interpreter", "latex", "fontsize", 10)
    xlabel("$M_{x}$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_{y}$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_{z}$", "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([0 1]);
    view(3)
    hold off
    % assign frame and write it to the video
    frame = getframe(h);
    video.writeVideo(frame);
end

video.close();