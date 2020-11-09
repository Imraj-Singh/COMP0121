clc
clear

%% Free procession with relaxation
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

% Define modelling parameters
RFPulse = 1;

% End time
t = 6;
resolution = 50;

% Duration and timestep of model
time = linspace(0, t, t * resolution + 1);

% The magnetic field strength caused by RF
B1 = dtheta/(gamma*RFPulse);

% Precession frequency
omega1 = gamma*B1;


% Forcing to 90 indices over 1 sec
endForce = ((length(time)-1)/6) + 1;


Msoln = zeros(length(time),3);

Msoln(1:endForce,1) = M(1);
Msoln(1:endForce,2) = M(2)*cos(omega1.*time(1:endForce)) + M(3)*sin(omega1.*time(1:endForce));
Msoln(1:endForce,3) = M(3)*cos(omega1.*time(1:endForce)) - M(2)*sin(omega1.*time(1:endForce));

endIndex = length(time);
timeRelax = time(endForce:endIndex)-1;
M = [0, 1, 0]';

Msoln(endForce:endIndex,1) = exp(-timeRelax./T2).*(M(1)*cos(omega0.*timeRelax) + M(2)*sin(omega0.*timeRelax));
Msoln(endForce:endIndex,2) = exp(-timeRelax./T2).*(M(2)*cos(omega0.*timeRelax) - M(1)*sin(omega0.*timeRelax));
Msoln(endForce:endIndex,3) = M(3)*exp(-timeRelax./T1) + 1*(1-exp(-timeRelax./T1));

%% Animation module

video = VideoWriter(['7_2', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = resolution;
video.set('FrameRate', frameRate);

video.open();

h = figure;

for i=1:length(time)
    quiver3(0,0,0,Msoln(i,1),Msoln(i,2),Msoln(i,3),'linewidth',5,'LineStyle','-','Color','r')
    hold on
    plot3(Msoln(1:i-1,1),Msoln(1:i-1,2),Msoln(1:i-1,3),'k--','linewidth',2)
    quiver3(0,0,0,0,0,1,'linewidth',2,'LineStyle','--','Color','k')
    title(['Time: ', num2str(time(i)),'s'], "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([0 1]);
    xlabel("$M_x$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_y$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_z$", "interpreter", "latex", "fontsize", 15)
    hold off
    frame = getframe(h);
    video.writeVideo(frame);
end

video.close();




