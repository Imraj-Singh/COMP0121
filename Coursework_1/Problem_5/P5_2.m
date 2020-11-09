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


% Calculate neccessary parameters

% Precession frequency
omega0 = 2 * 2 * pi;
deltaomega = 0;%omega0/100;
delta = .1;
omegaovector = omega0+deltaomega+delta*tan(pi*(rand(100,1)-1/2));

% Define modelling parameters
RFPulse = 1;

% The magnetic field strength caused by RF
B1 = dtheta/(gamma*RFPulse);

% Precession frequency
omega1 = gamma*B1;

dt = 0.1;

% Flip time 1
FT_1 = 0;

% Dephase time 1
DP_1 = 3;

% Flip time 2
FT_2 = 1;

% Rephase time
RP = 3;

% Echo time
Echo = 1;

% Dephase time 2
DP_2 = 3;

time = 0:dt:(FT_1 + DP_1 + FT_2 + RP + Echo + DP_2);

% Indices
% Flip time 1
FT_1i = FT_1/dt+1;

% Dephase time 1
DP_1i = (FT_1+DP_1)/dt+1;

% Flip time 2
FT_2i = (FT_1+DP_1+FT_2)/dt+1;

% Rephase time
RPi = (FT_1+DP_1+FT_2+RP)/dt+1;

% Echo time
Echoi = (FT_1+DP_1+FT_2+RP+Echo)/dt+1;

% Dephase time 2
DP_2i = length(time);


timeIll = 0:dt:1;
timeP = 0:dt:3;

Msoln = zeros(length(time),3, length(omegaovector));

for i = 1:length(omegaovector)
    % Flip time 1
    M = [0, 0, 1]';
    Msoln(1:FT_1i,1,i) = M(1);
    Msoln(1:FT_1i,2,i) = M(2);
    Msoln(1:FT_1i,3,i) = M(3);

    % Dephase time 1
    M = [0, 1, 0]';
    Msoln(FT_1i:DP_1i,1,i) = (M(1)*cos((omega0-omegaovector(i)).*timeP) + M(2)*sin((omega0-omegaovector(i)).*timeP));
    Msoln(FT_1i:DP_1i,2,i) = (M(2)*cos((omega0-omegaovector(i)).*timeP) - M(1)*sin((omega0-omegaovector(i)).*timeP));%.*exp(-timeRelax./T2);
    Msoln(FT_1i:DP_1i,3,i) = M(3);

    % Flip time 2
    M = Msoln(DP_1i,:,i);
    DP_1i = (FT_1+DP_1)/dt+1;
    Msoln(DP_1i:FT_2i,1,i) = M(1);
    Msoln(DP_1i:FT_2i,2,i) = M(2);
    Msoln(DP_1i:FT_2i,3,i) = M(3);

    % Rephase time
    Msoln(FT_2i:RPi,1,i) = (-M(1))*cos((omega0-omegaovector(i)).*timeP) + M(2)*sin((omega0-omegaovector(i)).*timeP);
    Msoln(FT_2i:RPi,2,i) = M(2)*cos((omega0-omegaovector(i)).*timeP) + M(1)*sin((omega0-omegaovector(i)).*timeP);%.*exp(-timeRelax./T2);
    Msoln(FT_2i:RPi,3,i) = M(3);

    % Echo time
    M = Msoln(RPi,:,i);
    Msoln(RPi:Echoi,1,i) = M(1);
    Msoln(RPi:Echoi,2,i) = M(2);%.*exp(-timeRelax./T2);
    Msoln(RPi:Echoi,3,i) = M(3);
    
    Echoi = (FT_1+DP_1+FT_2+RP+Echo)/dt+1;

    % Dephase time 2
    M = Msoln(Echoi,:,i);
    Msoln(Echoi:DP_2i,1,i) = M(1)*cos((omega0-omegaovector(i)).*timeP) + M(2)*sin((omega0-omegaovector(i)).*timeP);
    Msoln(Echoi:DP_2i,2,i) = M(2)*cos((omega0-omegaovector(i)).*timeP) - M(1)*sin((omega0-omegaovector(i)).*timeP);
    Msoln(Echoi:DP_2i,3,i) = M(3);
end
%% Animation module

video = VideoWriter(['5_2', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = 1/dt;
video.set('FrameRate', frameRate);

video.open();

h = figure;

for i=1:length(time)
    quiver3(0,0,0,Msoln(i,1,1),Msoln(i,2,1),Msoln(i,3,1),'linewidth',5,'LineStyle','-','Color','r')
    if time(i) > 0 && time(i) <= FT_1
        title('90$^\circ$ RF pulse around x-axis, instantaneous flip to +ve y-axis', "interpreter", "latex", "fontsize", 15)
    elseif time(i) > FT_1 && time(i) <= FT_1 + DP_1
        title(['Time: ', num2str(time(i) - FT_1),'s (Dephasing)'], "interpreter", "latex", "fontsize", 15)
    elseif time(i) > FT_1 + DP_1 && time(i) <= FT_1 + DP_1 + FT_2
        title('180$^\circ$ RF pulse around y-axis, $\tau = 3 s$', "interpreter", "latex", "fontsize", 15)
    elseif time(i) > FT_1 + DP_1 + FT_2 && time(i) <= FT_1 + DP_1 + FT_2 + RP
        title(['Time: ', num2str(time(i) - FT_1 - FT_2),'s (Rephasing)'], "interpreter", "latex", "fontsize", 15)
    elseif time(i) > FT_1 + DP_1 + FT_2 + RP && time(i) <= FT_1 + DP_1 + FT_2 + RP + Echo
        title('Echo time $T_E = 2 \tau = 6 s$', "interpreter", "latex", "fontsize", 15)
    elseif time(i) > FT_1 + DP_1 + FT_2 + RP + Echo && time(i) <= FT_1 + DP_1 + FT_2 + RP + Echo + DP_2
        title(['Time: ', num2str(time(i) - FT_1 - FT_2 - Echo),'s (Dephasing)'], "interpreter", "latex", "fontsize", 15)
    end
    %title(['Time: ', num2str(time(i)),'s'], "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([0 1]);
    xlabel("$M_x$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_y$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_z$", "interpreter", "latex", "fontsize", 15)
    hold on
    for j=2:length(omegaovector)
        quiver3(0,0,0,Msoln(i,1,j),Msoln(i,2,j),Msoln(i,3,j),'linewidth',5,'LineStyle','-','Color','r')
    end
    %plot3(Msoln(1:i-1,1),Msoln(1:i-1,2),Msoln(1:i-1,3),'b-')
    %quiver3(zeros(i-1,1),zeros(i-1,1),zeros(i-1,1),Msoln(1:i-1,1),Msoln(1:i-1,2),Msoln(1:i-1,3),'linewidth',1,'LineStyle','--','Color','b')
    hold off
    frame = getframe(h);
    video.writeVideo(frame);
end

video.close();




