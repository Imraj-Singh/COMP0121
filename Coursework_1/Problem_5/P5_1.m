clc
clear

%% Free procession with Dephasing
% Author: Imraj Singh 03/11/2020

% Given parameters

% Number of sample in each part of sequence
N = 101;

% Number of isochromats
Num = 1001;

% Prescribe T1
T2 = 1;
T1 = T2 * 2;

% Calculate neccessary parameters

% Precession frequency
omega0 = 50 * pi;
deltaomega = 16;
omegaovector = linspace(-deltaomega,deltaomega,Num);

% Inital time 1
INIT_1 = 0;
tINIT_1 = linspace(0, INIT_1, N);

% Flip time 1
FT_1 = .001;
tFT_1 = linspace(0, FT_1, N);

% Dephasing time 1
DP_1 = .5;
tDP_1 = linspace(FT_1, DP_1 + FT_1, N);

% Flip time 2
FT_2 = .001;
tFT_2 = linspace(DP_1 + FT_1, DP_1 + FT_1 + FT_2, N);

% Rephasing time
RP = DP_1;
tRP = linspace(DP_1 + FT_1 + FT_2, DP_1 + FT_1 + FT_2 + RP, N);

% Inital time 2
INIT_2 = 0;
tINIT_2 = linspace(DP_1 + FT_1 + FT_2 + RP,  DP_1 + FT_1 + FT_2 + RP, N);

% Dephasing time 2
DP_2 = DP_1;
tDP_2 = linspace(DP_1 + FT_1 + FT_2 + RP, + DP_1 + FT_1 + FT_2 + RP + DP_2, N);

time = [tINIT_1, tFT_1, tDP_1, tFT_2, tRP, tINIT_2, tDP_2];

% Indices
% Initial time
INITi_1 = N;

% Flip time 1
FTi_1 = INITi_1 + N;

% Dephasing time 1
DPi_1 = FTi_1 + N;

% Flip time 2
FTi_2 = DPi_1 + N;

% Rephasing time
RPi = FTi_2 + N;

% Initial time
INITi_2 = RPi + N;

% Dephasing time 2
DPi_2 = INITi_2 + N;


Msoln = zeros(length(time),3, length(omegaovector));

for i = 1:length(omegaovector)
    % Initial time
    Msoln(1:INITi_1,1,i) = 0;
    Msoln(1:INITi_1,2,i) = 0;
    Msoln(1:INITi_1,3,i) = 1;

    % Flip time 1
    tad = time;
    Vis90 = pi/2/FT_1;
    M = Msoln(INITi_1,:,i);
    Msoln(INITi_1:FTi_1,1,i) = M(1);
    Msoln(INITi_1:FTi_1,2,i) = M(2)*cos(Vis90.*tad(INITi_1:FTi_1)) + M(3)*sin(Vis90.*tad(INITi_1:FTi_1));
    Msoln(INITi_1:FTi_1,3,i) = M(3)*cos(Vis90.*tad(INITi_1:FTi_1)) - M(2)*sin(Vis90.*tad(INITi_1:FTi_1));

    % Dephase time 1 
    tad = tad - FT_1;
    M = Msoln(FTi_1,:,i);
    Msoln(FTi_1:DPi_1,1,i) = (M(1)*cos(omegaovector(i).*tad(FTi_1:DPi_1)) + M(2)*sin(omegaovector(i).*tad(FTi_1:DPi_1)));
    Msoln(FTi_1:DPi_1,2,i) = (M(2)*cos(omegaovector(i).*tad(FTi_1:DPi_1)) - M(1)*sin(omegaovector(i).*tad(FTi_1:DPi_1)));
    Msoln(FTi_1:DPi_1,3,i) = M(3);

    % Flip time 2
    tad = tad - DP_1;
    Vis90 = pi/FT_2;
    M = Msoln(DPi_1,:,i);
    Msoln(DPi_1:FTi_2,1,i) = M(1)*cos(Vis90.*tad(DPi_1:FTi_2)) + M(3)*sin(Vis90.*tad(DPi_1:FTi_2));
    Msoln(DPi_1:FTi_2,2,i) = M(2);
    Msoln(DPi_1:FTi_2,3,i) = M(3)*cos(Vis90.*tad(DPi_1:FTi_2)) - M(2)*sin(Vis90.*tad(DPi_1:FTi_2));

    % Rephase time
    tad = tad - FT_2;
    M = Msoln(FTi_2,:,i);
    Msoln(FTi_2:RPi,1,i) = M(1)*cos(omegaovector(i).*tad(FTi_2:RPi)) + M(2)*sin(omegaovector(i).*tad(FTi_2:RPi));
    Msoln(FTi_2:RPi,2,i) = M(2)*cos(omegaovector(i).*tad(FTi_2:RPi)) - M(1)*sin(omegaovector(i).*tad(FTi_2:RPi));
    Msoln(FTi_2:RPi,3,i) = M(3);

    % Echo time
    tad = tad - RP;
    M = Msoln(RPi,:,i);
    Msoln(RPi:INITi_2,1,i) = M(1);
    Msoln(RPi:INITi_2,2,i) = M(2);
    Msoln(RPi:INITi_2,3,i) = M(3);
    

    % Dephase time 2 
    tad = tad - INIT_2;
    M = Msoln(INITi_2,:,i);
    Msoln(INITi_2:DPi_2,1,i) = M(1)*cos(omegaovector(i).*tad(INITi_2:DPi_2)) + M(2)*sin(omegaovector(i).*tad(INITi_2:DPi_2));
    Msoln(INITi_2:DPi_2,2,i) = M(2)*cos(omegaovector(i).*tad(INITi_2:DPi_2)) - M(1)*sin(omegaovector(i).*tad(INITi_2:DPi_2));
    Msoln(INITi_2:DPi_2,3,i) = M(3);
end
Mperp = zeros(length(time),1);
for i=1:length(time)
    Mperp(i) = exp(-time(i)/T2)*sqrt(sum(Msoln(i,1,:))^2 + sum(Msoln(i,2,:))^2)/Num;
end

%% Animation module

video = VideoWriter(['5_1', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = N/8;
video.set('FrameRate', frameRate);

video.open();

h = figure;


for i=INITi_1:length(time)
    subplot(2,2,1)
    
    quiver3(0,0,0,Msoln(i,1,1),Msoln(i,2,1),Msoln(i,3,1),'linewidth',5,'LineStyle','-','Color','r')
    hold on
    for j=1:100:length(omegaovector)
        quiver3(0,0,0,Msoln(i,1,j),Msoln(i,2,j),Msoln(i,3,j),'linewidth',5,'LineStyle','-','Color','r')
    end
    grid on
    box on
    hold off
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([-1 1]);
    xlabel("$M_{x'}$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_{y'}$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_{z'}$", "interpreter", "latex", "fontsize", 15)
    
    subplot(2,2,2)
    quiver3(0,0,0,Msoln(i,1,1),Msoln(i,2,1),Msoln(i,3,1),'linewidth',5,'LineStyle','-','Color','r')
    hold on
    for j=1:100:length(omegaovector)
        quiver3(0,0,0,Msoln(i,1,j),Msoln(i,2,j),Msoln(i,3,j),'linewidth',5,'LineStyle','-','Color','r')
    end
    grid on
    box on
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([-1 1]);
    hold off
    view(2)
    xlabel("$M_{x'}$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_{y'}$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_{z'}$", "interpreter", "latex", "fontsize", 15)
    
    
    subplot(2,2,3:4)
    
    plot(time(1:i),Mperp(1:i),'linewidth',2,'LineStyle','-','Color','k')
    hold on
    plot(time(:),exp(-time(:)/T2),'linewidth',2,'LineStyle','--','Color','k')
    hold off
    grid on
    box on
    ylim([0 1]);
    xlim([0 max(time)]);
    xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    ylabel("Normalised signal", "interpreter", "latex", "fontsize", 15)
    legend("$\mid \frac{1}{N} \sum M_{\perp} \mid$","$T_2$ envelope", "interpreter", "latex", "fontsize", 10)
    
    % Initial
    if i > 0 && i <= INITi_1
        title('Magnetisation vector along $B_0$ initially', "interpreter", "latex", "fontsize", 15)
        
        % First flip
    elseif i > INITi_1 && i <= FTi_1
        xlim([0 FT_1]);
        title('RF pulse $\pi /2$ ', "interpreter", "latex", "fontsize", 15)
        
        % First Dephase
    elseif i > FTi_1 && i <= DPi_1
        title(['Isocromats dephasing, t = ', num2str(round(time(i) - INIT_1 - FT_1, 3)),' s'], "interpreter", "latex", "fontsize", 15)
        
        % Second flip
    elseif i > DPi_1 && i <= FTi_2
        title(['RF pulse $\pi$, $t = \tau$ = ', num2str(round(DP_1, 3)),' s'], "interpreter", "latex", "fontsize", 15)
    
        % Rephasing
    elseif i > FTi_2 && i <= RPi
        title(['Isocromats rephasing, $t =$ ', num2str(round(time(i) - INIT_1 - FT_1 - FT_2, 3)),' s (Relaxation)'], "interpreter", "latex", "fontsize", 15)
        
        % Echo
    elseif i > RPi && i <= INITi_2
        title(['Echo time, $t = T_E = 2 \tau =$ ',num2str(round(DP_1 * 2, 3)),' s'], "interpreter", "latex", "fontsize", 15)
        
        % Dephasing
    elseif i > INITi_2 && i <= DPi_2
        title(['Isocromats dephasing, $t =$ ', num2str(round(time(i) - INIT_1 - FT_1 - FT_2, 3)),' s (Relaxation)'], "interpreter", "latex", "fontsize", 15)
        
    end 
    
    frame = getframe(h);
    video.writeVideo(frame);
end

video.close();




