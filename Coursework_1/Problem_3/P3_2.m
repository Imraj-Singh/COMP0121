% housekeeping
clc
clear
%% Free procession ignoring relaxation
% Author: Imraj Singh 03/11/2020

% Given parameters
% Magnitisation aligned along y-axis initially
M = [0, 1, 0]';
M0 = [0, 0, 1]';

% Prescribe T1
T2 = 1;
T1 = T2 * 2;

% Calculate neccessary parameters

% Larmor frequency
omega0 = 2 * 2 * pi;

% End time
t = 6;
resolution = 50;

% Duration and timestep of model
time = linspace(0, t, t * resolution + 1);



Msoln(:,1) = exp(-time./T2).*(M(1)*cos(omega0.*time) + M(2)*sin(omega0.*time));
Msoln(:,2) = exp(-time./T2).*(M(2)*cos(omega0.*time) - M(1)*sin(omega0.*time));
Msoln(:,3) = M(3)*exp(-time./T1) + M0(3)*(1-exp(-time./T1));


%% Animation module
% Author: Imraj Singh 03/11/2020

% initialise the video
video = VideoWriter(['3_2', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = resolution;
video.set('FrameRate', frameRate);

video.open();

% initialise the figure
h = figure;

% specify animation captures each degree of rotation
for i=1:length(time)
    quiver3(0,0,0,Msoln(i,1),Msoln(i,2),Msoln(i,3),'linewidth',4,'LineStyle','-','Color','r')
    hold on
    plot3(Msoln(1:i-1,1),Msoln(1:i-1,2),Msoln(1:i-1,3),'k--','linewidth',2)
    quiver3(0,0,0,0,1,0,'linewidth',2,'LineStyle','--','Color','k')

    % format title, legend, labels, limits, grid and box
    title(['Time: ', num2str(time(i)),' s'], "interpreter", "latex", "fontsize", 15)
    xlabel("$M_x$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_y$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_z$", "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([0 1]);
    hold off
    
    % assign frame and write it to the video
    frame = getframe(h);
    video.writeVideo(frame);
end

video.close();

