clc
clear

%% Problem 4.1 IR sequence magnetisation vector visualisation
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

% Flip time 1
FT_1 = .002;
tFT_1 = linspace(0, FT_1, N);

% Relaxation time 1
RT_1 = .5;
tRT_1 = linspace(FT_1, RT_1 + FT_1, N);

% Flip time 2
FT_2 = .001;
tFT_2 = linspace(RT_1 + FT_1, RT_1 + FT_1 + FT_2, N);

% Relaxation time 2
RT_2 = 8;
tRT_2 = linspace(RT_1 + FT_1 + FT_2, + RT_1 + FT_1 + FT_2 + RT_2, N);

time = [tINIT, tFT_1, tRT_1, tFT_2, tRT_2];

% Indices
% Initial time
INITi_1 = N;

% Flip time 1
FTi_1 = INITi_1 + N;

% Relaxation time 1
RTi_1 = FTi_1 + N;

% Flip time 2
FTi_2 = RTi_1 + N;

% Relaxation time 2
RTi_2 = FTi_2 + N;

Msoln = zeros(length(time),3);

% Initial time
Msoln(1:INITi_1,1) = M(1);
Msoln(1:INITi_1,2) = M(2);
Msoln(1:INITi_1,3) = M(3);

% Flip time 1
tad = time;
Vis90 = pi/FT_1;
M = Msoln(INITi_1,:);
Msoln(INITi_1:FTi_1,1) = M(1);
Msoln(INITi_1:FTi_1,2) = M(2)*cos(Vis90.*tad(INITi_1:FTi_1)) + M(3)*sin(Vis90.*tad(INITi_1:FTi_1));
Msoln(INITi_1:FTi_1,3) = M(3)*cos(Vis90.*tad(INITi_1:FTi_1)) - M(2)*sin(Vis90.*tad(INITi_1:FTi_1));


% Relaxation time 1
tad = tad - FT_1;
M = Msoln(FTi_1,:);
Msoln(FTi_1:RTi_1,1) = M(1).*exp(-tad(FTi_1:RTi_1)./T2);
Msoln(FTi_1:RTi_1,2) = M(2).*exp(-tad(FTi_1:RTi_1)./T2);
Msoln(FTi_1:RTi_1,3) = M(3)*exp(-tad(FTi_1:RTi_1)./T1) + (1-exp(-tad(FTi_1:RTi_1)./T1));

% Flip time 2
tad = tad - RT_1;
Vis90 = pi/2/FT_2;
M = Msoln(RTi_1,:);
Msoln(RTi_1:FTi_2,1) = M(1);
Msoln(RTi_1:FTi_2,2) = M(2)*cos(Vis90.*tad(RTi_1:FTi_2)) + M(3)*sin(Vis90.*tad(RTi_1:FTi_2));
Msoln(RTi_1:FTi_2,3) = M(3)*cos(Vis90.*tad(RTi_1:FTi_2)) - M(2)*sin(Vis90.*tad(RTi_1:FTi_2));

% Relaxation time 2
tad = tad - FT_2;
M = Msoln(FTi_2,:);
Msoln(FTi_2:RTi_2,1) = M(1).*exp(-tad(FTi_2:RTi_2)./T2);
Msoln(FTi_2:RTi_2,2) = M(2).*exp(-tad(FTi_2:RTi_2)./T2);
Msoln(FTi_2:RTi_2,3) = M(3)*exp(-tad(FTi_2:RTi_2)./T1) + (1-exp(-tad(FTi_2:RTi_2)./T1));
    

%% Animation module
% Author: Imraj Singh 03/11/2020

% initialise the video
video = VideoWriter(['4_2_rot_opt', '.mp4'], 'MPEG-4');

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
    zlim([-1 1]);
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
    ylim([-1 1]);
    hold off
    
    subplot(2,2,3:4)
    plot(time(1:i),Msoln(1:i,1),'linewidth',2,'LineStyle','-','Color','r');
    hold on
    plot(time(1:i),Msoln(1:i,2),'linewidth',2,'LineStyle','-.','Color','b');
    plot(time(1:i),Msoln(1:i,3),'linewidth',2,'LineStyle','--','Color','k');
    
    % format title, legend, labels, limits, grid and box
    
    if i > 0 && i <= INITi_1
        title('Magnetisation vector along $B_0$ initially', "interpreter", "latex", "fontsize", 15)
        xlim([0 1]);
        xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    elseif i > INITi_1 && i <= FTi_1
        xlim([0 FT_1]);
        title(['RF pulse $\pi $, time: ',num2str(round(1000*time(i), 3)),' ms'], "interpreter", "latex", "fontsize", 15)
        xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    elseif i > FTi_1 && i <= RTi_1
        title(['Time from $\pi$ flip: ', num2str(round(time(i) - INIT - FT_1, 3)),' s (Relaxation)'], "interpreter", "latex", "fontsize", 15)
        xlim([0 FT_1+RT_1+FT_2+RT_2]);
        xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    elseif i > RTi_1 && i <= FTi_2
        title(['RF pulse $\pi / 2$, time: ',num2str(round(1000*(time(i) - FT_1 - RT_1), 3)),' ms'], "interpreter", "latex", "fontsize", 15)
        xlim([0 FT_1+RT_1+FT_2+RT_2]);
        xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    elseif i > FTi_2 && i <= RTi_2
        title(['Time from $\pi / 2$ flip: ', num2str(round(time(i) - INIT - FT_1 - RT_1 - FT_2, 3)),' s (Relaxation)'], "interpreter", "latex", "fontsize", 15)
        xlim([0 FT_1+RT_1+FT_2+RT_2]);
        xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    end
    
    ylabel("Magnetisation", "interpreter", "latex", "fontsize", 15)
    legend("$M_{x'}$","$M_{y'}$", "$M_{z'}$","interpreter", "latex", "fontsize", 10)
    grid on
    box on
    ylim([-1 1]);
    hold off
    
    
    % assign frame and write it to the video
    frame = getframe(h);
    video.writeVideo(frame);
end

video.close();


