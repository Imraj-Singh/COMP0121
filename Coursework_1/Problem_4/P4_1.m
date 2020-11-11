clc
clear

%% Problem 4.1 FID sequence magnetisation vector visualisation
% Author: Imraj Singh 11/11/2020

% Given parameters

% Magnitisation aligned along z-axis initially
M = [0, 0, 1]';

% Number of sample in each part of sequence
N = 101;

% Prescribe T1
T2 = 1;
T1 = T2 * 2;

% Calculate neccessary parameters

% Precession frequency
omega0 = 2 * 2 * pi;

% Flip time
INIT = 0;
tINIT = linspace(0, 0, N);

% Flip time
FT = .001;
tFT = linspace(0, FT, N);

% Relaxation time
RT = 8;
tRT = linspace(FT, RT + FT, N);

time = [tINIT, tFT, tRT];

% Indices
% Initial time
INITi = N;

% Flip time
FTi = INITi + N;

% Relaxation time
RTi = FTi + N;

Msoln = zeros(length(time),3);

% Initial time
Msoln(1:INITi,1) = M(1);
Msoln(1:INITi,2) = M(2);
Msoln(1:INITi,3) = M(3);

% Flip time
Vis90 = pi/2/FT;
M = Msoln(INITi,:);
Msoln(INITi:FTi,1) = M(1);
Msoln(INITi:FTi,2) = M(2)*cos(Vis90.*time(INITi:FTi)) + M(3)*sin(Vis90.*time(INITi:FTi));
Msoln(INITi:FTi,3) = M(3)*cos(Vis90.*time(INITi:FTi)) - M(2)*sin(Vis90.*time(INITi:FTi));

% Flip time 2
M = Msoln(FTi,:);
Msoln(FTi:RTi,1) = M(1).*exp(-time(FTi:RTi)./T2);
Msoln(FTi:RTi,2) = M(2).*exp(-time(FTi:RTi)./T2);
Msoln(FTi:RTi,3) = M(3)*exp(-time(FTi:RTi)./T1) + (1-exp(-time(FTi:RTi)./T1));
    

%% Animation module
% Author: Imraj Singh 03/11/2020

% initialise the video
video = VideoWriter(['4_1_rot', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = (N-1)/10;
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

    % format legend, labels, limits, grid and box
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
    quiver(0,0,Msoln(i,2),Msoln(i,3),'linewidth',4,'LineStyle','-','Color','r')
    hold on
    plot(Msoln(1:i,2),Msoln(1:i,3),'k--','linewidth',2)
    quiver(0,0,0,1,'linewidth',2,'LineStyle','--','Color','k')
    
    % format legend, labels, limits, grid and box
    xlabel("$M_{y'}$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_{z'}$", "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([0 1]);
    hold off
    
    subplot(2,2,3:4)
    plot(time(1:i),Msoln(1:i,1),'linewidth',2,'LineStyle','-','Color','r');
    hold on
    plot(time(1:i),Msoln(1:i,2),'linewidth',2,'LineStyle','-.','Color','b');
    plot(time(1:i),Msoln(1:i,3),'linewidth',2,'LineStyle','--','Color','k');
    
    % format title, legend, labels, limits, grid and box
    
    if i > 0 && i <= INITi
        title('Magnetisation vector along $B_0$ initially', "interpreter", "latex", "fontsize", 15)
    xlim([0 1]);
    xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    elseif i > INITi && i <= FTi
        title(['RF pulse $\pi / 2$, time: ',num2str(round(1000*time(i), 3)),' ms'], "interpreter", "latex", "fontsize", 15)
    xlim([0 FT]);
    xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    elseif i > FTi && i <= RTi
        title(['Time from $\pi / 2$ flip: ', num2str(round(time(i) - INIT - FT, 3)),' s (Relaxation)'], "interpreter", "latex", "fontsize", 15)
    xlim([0 RT]);
    xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    end
    ylabel("Magnetisation", "interpreter", "latex", "fontsize", 15)
    legend("$M_{x'}$","$M_{y'}$", "$M_{z'}$","interpreter", "latex", "fontsize", 10)
    grid on
    box on
    ylim([0 1]);
    hold off
    
    
    % assign frame and write it to the video
    frame = getframe(h);
    video.writeVideo(frame);
end

video.close();


