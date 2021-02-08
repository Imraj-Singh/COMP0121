% housekeeping
clc
clear
%% Free procession with relaxation rotating FoR
% Author: Imraj Singh 03/11/2020

% Given parameters
% Magnitisation aligned along y-axis initially
M = [0, 1, 0]';
M0 = [0, 0, 1]';

% Prescribe T1
T2 = 1;
T1 = T2 * 2;

% Calculate neccessary parameters

% Precession frequency
omega0 = 2 * 2 * pi;

% End time
t = 6;
resolution = 50;

% Duration and timestep of model
time = linspace(0, t, t * resolution + 1);


Msoln(:,1) = exp(-time./T2).*M(1);
Msoln(:,2) = exp(-time./T2).*M(2);
Msoln(:,3) = M(3)*exp(-time./T1) + M0(3)*(1-exp(-time./T1));


%% Animation

    subplot(2,2,1)
    quiver3(0,0,0,0,1,0,'linewidth',2,'LineStyle','--','Color','k')
    hold on
    plot3(Msoln(:,1),Msoln(:,2),Msoln(:,3),'k--','linewidth',2)
    quiver3(0,0,0,Msoln(end,1),Msoln(end,2),Msoln(end,3),'linewidth',4,'LineStyle','-','Color','r')

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
    quiver(0,0,Msoln(end,2),Msoln(end,3),'linewidth',4,'LineStyle','-','Color','r')
    hold on
    plot(Msoln(:,2),Msoln(:,3),'k--','linewidth',2)
    quiver(0,0,1,0,'linewidth',2,'LineStyle','--','Color','k')
    
    % format legend, labels, limits, grid and box
    xlabel("$M_{y'}$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_{z'}$", "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([0 1]);
    hold off
    
    subplot(2,2,3:4)
    plot(time(:),Msoln(:,1),'linewidth',2,'LineStyle','-','Color','r');
    hold on
    plot(time(:),Msoln(:,2),'linewidth',2,'LineStyle','-.','Color','b');
    plot(time(:),Msoln(:,3),'linewidth',2,'LineStyle','--','Color','k');
    
    % format title, legend, labels, limits, grid and box
    xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    ylabel("Magnetisation", "interpreter", "latex", "fontsize", 15)
    legend("$M_{x'}$","$M_{y'}$", "$M_{z'}$","interpreter", "latex", "fontsize", 10)
    grid on
    box on
    xlim([0 max(time)]);
    ylim([0 1]);
    hold off
    
    


