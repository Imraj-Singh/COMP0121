clc
clear

%% Problem 4.1 FID sequence magnetisation vector visualisation
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

dt = 0.05;

% Flip time
INIT = 1;

% Flip time
FT = 1;

% Relaxation time
RT = 8;

time = 0:dt:(INIT + FT + RT);

tINIT = 0:dt:INIT;
tFT = 0:dt:FT;
tRT = 0:dt:RT;

% Indices
% Initial time
INITi = INIT/dt+1;

% Flip time
FTi = (INIT + FT)/dt+1;

% Relaxation time
RTi = (INIT + FT + RT)/dt+1;

Msoln = zeros(length(time),3);

% Initial time
M = [0, 0, 1]';
Msoln(1:INITi,1) = M(1);
Msoln(1:INITi,2) = M(2);
Msoln(1:INITi,3) = M(3);

% Flip time
Vis90 = pi/2/FT;
M = Msoln(INITi,:);
Msoln(INITi:FTi,1) = M(1);
Msoln(INITi:FTi,2) = M(2)*cos(Vis90.*tFT) + M(3)*sin(Vis90.*tFT);
Msoln(INITi:FTi,3) = M(3)*cos(Vis90.*tFT) - M(2)*sin(Vis90.*tFT);

% Flip time 2
M = Msoln(FTi,:);
Msoln(FTi:RTi,1) = M(1).*exp(-tRT./T2);
Msoln(FTi:RTi,2) = M(2).*exp(-tRT./T2);
Msoln(FTi:RTi,3) = M(3)*exp(-tRT./T1) + (1-exp(-tRT./T1));

%% Animation module

video = VideoWriter(['4_1rot', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = 1/dt/4;
video.set('FrameRate', frameRate);

video.open();

h = figure;

for i=1:length(time)
    if time(i) > 0 && time(i) <= INIT
        quiver(0,0,Msoln(i,2),Msoln(i,3),'linewidth',5,'LineStyle','-','Color','r')
        hold on
        title('Magnetisation vector along $B_0$ initially', "interpreter", "latex", "fontsize", 15)
    elseif time(i) > INIT && time(i) <= INIT + FT
        quiver(0,0,Msoln(INITi,2),Msoln(INITi,3),'linewidth',2,'LineStyle','--','Color','k')
        hold on
        quiver(0,0,Msoln(FTi,2),Msoln(FTi,3),'linewidth',5,'LineStyle','-','Color','r')
        plot(Msoln(INITi:FTi,2),Msoln(INITi:FTi,3),'k--','linewidth',2)
        title('Instantaneous $\pi / 2$ flip', "interpreter", "latex", "fontsize", 15)
    elseif time(i) > INIT + FT && time(i) <= INIT + FT + RT
        plot(Msoln(1:i-1,2),Msoln(1:i-1,3),'k--','linewidth',2)
        hold on
        quiver(0,0,0,1,'linewidth',2,'LineStyle','--','Color','k')
        quiver(0,0,Msoln(i,2),Msoln(i,3),'linewidth',5,'LineStyle','-','Color','r')
        title(['Time from $\pi / 2$ flip: ', num2str(time(i) - INIT - FT),'s (Relaxation)'], "interpreter", "latex", "fontsize", 15)
    end
    grid on
    box on
    xlim([-1 1]);
    ylim([-1 1]);
    xlabel("$M_{y'}$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_{z'}$", "interpreter", "latex", "fontsize", 15)
    hold off
    frame = getframe(h);
    video.writeVideo(frame);
end

video.close();




