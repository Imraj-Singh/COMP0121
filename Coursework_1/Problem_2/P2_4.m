% housekeeping
clc
clear
%% Forced procession ignoring relaxation
% Author: Imraj Singh 03/11/2020

% Given parameters

% Magnitisation aligned along z-axis initially
M = [0, 0, 1]';

% Duration of the RF pulse
t = 1*10^(-3);

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

%% Animation

% initialise the video
video = VideoWriter(['2_4', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = 10;
video.set('FrameRate', frameRate);

video.open();

% initialise the figure
h = figure;

% specify animation captures each degree of rotation
for i=1:length(time)
    subplot(2,2,1)
    quiver3(0,0,0,0,0,1,'linewidth',2,'LineStyle','--','Color','k')
    hold on
    plot3(Msoln(1:i,1),Msoln(1:i,2),Msoln(1:i,3),'k--','linewidth',2)
    quiver3(0,0,0,Msoln(i,1),Msoln(i,2),Msoln(i,3),'linewidth',4,'LineStyle','-','Color','r')
    
    % format title, legend, labels, limits, grid and box
    xlabel("$M_{x'}$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_{y'}$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_{z'}$", "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([0 1]);
    view(3)
    hold off
    
    
    subplot(2,2,2)
    quiver(0,0,0,1,'linewidth',2,'LineStyle','--','Color','k')
    hold on
    plot(Msoln(1:i,2),Msoln(1:i,3),'k--','linewidth',2)
    quiver(0,0,Msoln(i,2),Msoln(i,3),'linewidth',4,'LineStyle','-','Color','r')
    
    % format title, legend, labels, limits, grid and box
    xlabel("$M_{y'}$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_{z'}$", "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([0 1]);
    hold off
    
    subplot(2,2,3:4)
    plot(time(1:i)*1000,Msoln(1:i,2),'linewidth',2,'LineStyle','-','Color','r');
    hold on
    plot(time(1:i)*1000,Msoln(1:i,3),'linewidth',2,'LineStyle','-.','Color','b');
    
    % format title, legend, labels, limits, grid and box
    title(['Time: ', num2str(time(i)*1000),' ms'], "interpreter", "latex", "fontsize", 15)
    xlabel("Time (ms)", "interpreter", "latex", "fontsize", 15)
    ylabel("Magnetisation", "interpreter", "latex", "fontsize", 15)
    legend("$M_{y'}$","$M_{z'}$", "interpreter", "latex", "fontsize", 10)
    grid on
    box on
    xlim([0 max(time*1000)]);
    ylim([0 1]);
    hold off
    
    % assign frame and write it to the video
    frame = getframe(h);
    video.writeVideo(frame);
end

video.close();



