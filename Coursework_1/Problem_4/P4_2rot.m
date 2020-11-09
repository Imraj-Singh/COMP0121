clc
clear

%% Problem 4.1 FID sequence magnetisation vector visualisation
% Author: Imraj Singh 03/11/2020

% Given parameters

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
FT_1 = 1;

% Relaxation time
RT_1 = 0.5;

% Flip time
FT_2 = 1;

% Relaxation time
RT_2 = 8;

time = 0:dt:(INIT + FT_1 + RT_1 + FT_2 + RT_2);

tINIT = 0:dt:INIT;
tFT_1 = 0:dt:FT_1;
tRT_1 = 0:dt:RT_1;
tFT_2 = 0:dt:FT_2;
tRT_2 = 0:dt:RT_2;

% Indices
% Initial time
INITi = INIT/dt+1;

% Flip time
FTi_1 = (INIT + FT_1)/dt+1;

% Relaxation time
RTi_1 = (INIT + FT_1 + RT_1)/dt+1;

% Flip time
FTi_2 = (INIT + FT_1 + RT_1 + FT_2)/dt+1;

% Relaxation time
RTi_2 = (INIT + FT_1 + RT_1 + FT_2 + RT_2)/dt+1;

Msoln = zeros(length(time),3);

% Initial time
M = [0, 0, 1]';
Msoln(1:INITi,1) = M(1);
Msoln(1:INITi,2) = M(2);
Msoln(1:INITi,3) = M(3);

% Flip time 1
Vis90 = pi/FT_1;
M = Msoln(INITi,:);
Msoln(INITi:FTi_1,1) = M(1);
Msoln(INITi:FTi_1,2) = M(2)*cos(Vis90.*tFT_1) + M(3)*sin(Vis90.*tFT_1);
Msoln(INITi:FTi_1,3) = M(3)*cos(Vis90.*tFT_1) - M(2)*sin(Vis90.*tFT_1);

% Relaxation time 1
M = Msoln(FTi_1,:);
Msoln(FTi_1:RTi_1,1) = M(1).*exp(-tRT_1./T2);
Msoln(FTi_1:RTi_1,2) = M(2).*exp(-tRT_1./T2);
Msoln(FTi_1:RTi_1,3) = M(3)*exp(-tRT_1./T1) + (1-exp(-tRT_1./T1));

% Flip time 2
Vis90 = pi/FT_2/2;
M = Msoln(RTi_1,:);
Msoln(RTi_1:FTi_2,1) = M(1);
Msoln(RTi_1:FTi_2,2) = M(2)*cos(Vis90.*tFT_2) + M(3)*sin(Vis90.*tFT_2);
Msoln(RTi_1:FTi_2,3) = M(3)*cos(Vis90.*tFT_2) - M(2)*sin(Vis90.*tFT_2);

% Relaxation time 2
M = Msoln(FTi_2,:);
Msoln(FTi_2:RTi_2,1) = M(1).*exp(-tRT_2./T2);
Msoln(FTi_2:RTi_2,2) = M(2).*exp(-tRT_2./T2);
Msoln(FTi_2:RTi_2,3) = M(3)*exp(-tRT_2./T1) + (1-exp(-tRT_2./T1));

%% Animation module

video = VideoWriter(['4_2rot', '.mp4'], 'MPEG-4');

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
    elseif time(i) > INIT && time(i) <= INIT + FT_1
        quiver(0,0,Msoln(INITi,2),Msoln(INITi,3),'linewidth',2,'LineStyle','--','Color','k')
        hold on
        quiver(0,0,Msoln(FTi_1,2),Msoln(FTi_1,3),'linewidth',5,'LineStyle','-','Color','r')
        plot(Msoln(INITi:FTi_1,2),Msoln(INITi:FTi_1,3),'k--','linewidth',2)
        title('Instantaneous $\pi$ flip', "interpreter", "latex", "fontsize", 15)
    elseif time(i) > INIT + FT_1 && time(i) <= INIT + FT_1 + RT_1
        plot(Msoln(1:i-1,2),Msoln(1:i-1,3),'k--','linewidth',2)
        hold on
        quiver(0,0,0,1,'linewidth',2,'LineStyle','--','Color','k')
        quiver(0,0,Msoln(i,2),Msoln(i,3),'linewidth',5,'LineStyle','-','Color','r')
        title(['Time from $\pi$ flip: ', num2str(time(i) - INIT - FT_1),'s (Relaxation)'], "interpreter", "latex", "fontsize", 15)
    elseif time(i) > INIT + FT_1 + RT_1 && time(i) <= INIT + FT_1 + RT_1 + FT_2
        plot(Msoln(1:FTi_2,2),Msoln(1:FTi_2,3),'k--','linewidth',2)
        hold on
        quiver(0,0,0,1,'linewidth',2,'LineStyle','--','Color','k')
        quiver(0,0,Msoln(FTi_2,2),Msoln(FTi_2,3),'linewidth',5,'LineStyle','-','Color','r')
        title('Instantaneous $\pi / 2$ flip', "interpreter", "latex", "fontsize", 15)
    elseif time(i) > INIT + FT_1 + RT_1 + FT_2 && time(i) <= INIT + FT_1 + RT_1 + FT_2 + RT_2
        plot(Msoln(1:i-1,2),Msoln(1:i-1,3),'k--','linewidth',2)
        hold on
        quiver(0,0,0,1,'linewidth',2,'LineStyle','--','Color','k')
        quiver(0,0,Msoln(i,2),Msoln(i,3),'linewidth',5,'LineStyle','-','Color','r')
        title(['Time from $\pi / 2$ flip: ', num2str(time(i) - INIT - FT_1 - FT_2),'s (Relaxation)'], "interpreter", "latex", "fontsize", 15)
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




