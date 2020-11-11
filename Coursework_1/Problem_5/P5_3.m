% house keeping
clc
clear
%% Sensitivity of Uniform Distribution to delta omega
% Author: Imraj Singh


% Given parameters

% Number of sample in each part of sequence
N = 101;

% Number of isochromats
Num = 10001;

% Prescribe T1
T2 = 1;
T1 = T2 * 2;

% Calculate neccessary parameters

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

%% Sensitivity of Uniform Distribution to delta omega
% Author: Imraj Singh

% Precession frequency
omega0 = 50 * pi;
delta = [1 2 4 8 16 32 64 128];

Mperp = zeros(length(time),length(delta));
for Di = 1:length(delta)
    omegaovector = linspace(-delta(Di),delta(Di),Num);
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
        Msoln(DPi_1:FTi_2,1,i) = M(1);
        Msoln(DPi_1:FTi_2,2,i) = M(2)*cos(Vis90.*tad(DPi_1:FTi_2)) + M(3)*sin(Vis90.*tad(DPi_1:FTi_2));
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
    for i=1:length(time)
        Mperp(i,Di) = exp(-time(i)/T2)*sqrt(sum(Msoln(i,1,:))^2 + sum(Msoln(i,2,:))^2)/Num;
    end
end

video = VideoWriter(['5_3_dOmega', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = .25;
video.set('FrameRate', frameRate);

video.open();

h = figure;
for i=1:length(delta)
    plot(time(:),Mperp(:,i),'linewidth',2,'LineStyle','-','Color','k')
    
    grid on
    box on
    xlim([0 max(time)]);
    ylim([0 1]);
    view(2)
    title(['$\delta \omega $ = ', num2str(delta(i))], "interpreter", "latex", "fontsize", 15)
    xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    ylabel("$\mid \frac{1}{N} \sum M_{\perp} \mid$", "interpreter", "latex", "fontsize", 15)
    
    frame = getframe(h);
    video.writeVideo(frame);
end

video.close();

%% Sensitivity of Cauchy to delta
% Author: Imraj Singh

% Precession frequency
omega0 = 50 * pi;
delta = [1 2 4 8 16 32 64 128];
deltaomega = 0;
Mperp = zeros(length(time),length(delta));
for Di = 1:length(delta)
    omegaovector = deltaomega+delta(Di)*tan(pi*(rand(Num,1)-1/2));
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
        Msoln(DPi_1:FTi_2,1,i) = M(1);
        Msoln(DPi_1:FTi_2,2,i) = M(2)*cos(Vis90.*tad(DPi_1:FTi_2)) + M(3)*sin(Vis90.*tad(DPi_1:FTi_2));
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
    for i=1:length(time)
        Mperp(i,Di) = exp(-time(i)/T2)*sqrt(sum(Msoln(i,1,:))^2 + sum(Msoln(i,2,:))^2)/Num;
    end
end

h = figure;
for i=1:length(delta)
    plot(time(:),Mperp(:,i),'linewidth',2,'LineStyle','-','Color','k')
    
    grid on
    box on
    xlim([0 max(time)]);
    ylim([0 1]);
    view(2)
    title(['$\Delta $ = ', num2str(delta(i))], "interpreter", "latex", "fontsize", 15)
    xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    ylabel("$\mid \frac{1}{N} \sum M_{\perp} \mid$", "interpreter", "latex", "fontsize", 15)
    
    frame(i) = getframe(h);
end
video = VideoWriter(['5_3_Delta', '.mp4'], 'MPEG-4');
frameRate = .25;
video.set('FrameRate', frameRate);
video.open();
video.writeVideo(frame);
video.close();

%% Sensitivity of Cauchy to omega
% Author: Imraj Singh

% Precession frequency
omega0 = 50 * pi;
delta1 = 4;
deltaomega = [-100 -50 -25 -10 -5 0 5 10 25 50 100];
Mperp = zeros(length(time),length(delta));
for Di = 1:length(delta)
    omegaovector = delta(Di)+delta1*tan(pi*(rand(Num,1)-1/2));
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
        Msoln(DPi_1:FTi_2,1,i) = M(1);
        Msoln(DPi_1:FTi_2,2,i) = M(2)*cos(Vis90.*tad(DPi_1:FTi_2)) + M(3)*sin(Vis90.*tad(DPi_1:FTi_2));
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
    for i=1:length(time)
        Mperp(i,Di) = exp(-time(i)/T2)*sqrt(sum(Msoln(i,1,:))^2 + sum(Msoln(i,2,:))^2)/Num;
    end
end


% set the frame rate

h = figure;
for i=1:length(delta)
    plot(time(:),Mperp(:,i),'linewidth',2,'LineStyle','-','Color','k')
    
    grid on
    box on
    xlim([0 max(time)]);
    ylim([0 1]);
    view(2)
    title(['$\omega_0 $ = ', num2str(delta(i))], "interpreter", "latex", "fontsize", 15)
    xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    ylabel("$\mid \frac{1}{N} \sum M_{\perp} \mid$", "interpreter", "latex", "fontsize", 15)
    frame(i) = getframe(gcf);
end
video = VideoWriter(['5_3_omega', '.mp4'], 'MPEG-4');
frameRate = .25;
video.set('FrameRate', frameRate);
video.open();
video.writeVideo(frame);
video.close();