%% Free procession with Dephasing
% Author: Imraj Singh 03/11/2020

clc
clear
close all

% Add the path what holds functions
addpath(genpath('.../Functions'))

%% Define coordinate position of each isochromat

% Number of isochromats x
Numx = 21;

% Number of isochromats y
Numy = 1;

% Spatial coordinates of isochromats (m)
x = linspace(-10,10,Numx)/1000;
y = linspace(0,0,Numy)/1000;

% Create mesh grid of x and y coordinates
[X, Y] = meshgrid(x,y);

% Flatten the array into a vector to make manipulation easier
X = X(:);
Y = Y(:);

%% Set NMR parameters and initialise isochromat structures

% Gyromagnetic ratio of hyrodrogen proton
gamma = 2.68*10^8;

% Initalise the isochromats, setting the X, Y cooridates and instantiating
% the Mx, My, Mz fields of the structures
iso = initalise_iso(X, Y);

% B0
B0 = 3;

% Initalise the B0 this is done such that arbitrary perturbation can be
% added to the B0 field to simulate field inhomogenities
iso = initalise_B0(iso, B0);

% Larmor frequency of rotating frame of reference
R_FoR = B0*gamma;

% Setting the brain tissure relaxation values. Approximation based off:
% Wansapura JP, Holland SK, Dunn RS, Ball WS Jr. NMR relaxation times in
% the human brain at 3.0 tesla. J Magn Reson Imaging. 1999 Apr;9(4):531-8.
% doi: 10.1002/(sici)1522-2586(199904)9:4<531::aid-jmri4>3.0.co;2-l.
% PMID: 10232510. Approximate average of white and grey matter.
T1 = 1;
T2 = .1;

%% Define the experiment in terms of blocks

% Each experiment begins with initialising the direction of the spins then
% by flipping the spins with a pi/2 hard RF pulse along the x' axis over
% 0.32ms

% [Set the number of points to be calculated minimum of 2 for each block
blocks.N = [3 11 11 101];

% Time over which the block occurs
blocks.t = [0 .32 5.12 512]/1000;

% The gradient of the field dB/dx for each block
blocks.Gx = [0 0 0 0]/1000;

% The gradient of the field dB/dz for each block
blocks.Gy = [0 0 0 0]/1000;

% The relaxation time only when spins are relaxation
blocks.relax = [0 0 5.12 512]/1000;

% Initialise the Mx, My, Mz values that stores the sum of all the spins
% magnetisation
Mx = zeros(sum(blocks.N),1);
My = zeros(sum(blocks.N),1);
Mz = zeros(sum(blocks.N),1);

% Calculate magnetisation for all the spins
for i=1:length(iso)
    % Initialise spin
    iso(i) = initial_dir(iso(i),'z',blocks.N(1));
    
    % RF pulse
    iso(i) = rf_pulse(iso(i), 'x', pi/2, blocks.N(2), blocks.t(2));
    
    % Free relaxation no gradient
    iso(i) = free_relaxation(iso(i), blocks.N(3), blocks.t(3), ...
        blocks.Gx(3), blocks.Gy(3), blocks.relax(1:3), gamma, T1, T2, R_FoR);
    
    % Free relaxation no gradient
    iso(i) = free_relaxation(iso(i), blocks.N(4), blocks.t(4), ...
        blocks.Gx(4), blocks.Gy(4), blocks.relax(1:4), gamma, T1, T2, R_FoR);
    
    % Add the magnetisations to the total  magnetisation
    Mx = Mx + iso(i).Mx;
    My = My + iso(i).My;
    Mz = Mz + iso(i).Mz;
end

% Magnetisation now normalised across all the spins
Mx = Mx/length(iso);
My = My/length(iso);
Mz = Mz/length(iso);

% Calculate the time vector
Time = unroll_time(blocks.N,blocks.t);

% Calculate the Gx, Gy, kx, ky vectors in time
[Gx_t, Gy_t, kx, ky] = unroll_plotting(blocks.N,blocks.Gx,blocks.Gy,gamma,Time);

%% Animation module

norm = 0.001;
norm_scale = 2;

cmap = colormap(hsv(360));
for i=1:length(Time)
    subplot(2,2,1)
    h = quiver3(iso(1).x,iso(1).y,0,iso(1).Mx(i)*norm,iso(1).My(i)*norm,iso(1).Mz(i)*norm,'linewidth',2,'LineStyle','-');
    alpha = atan2(iso(1).My(i)*norm, iso(1).Mx(i)*norm);
    idx = ceil(rad2deg(alpha));
    if idx < 1 % negative indices and 0 are invalid
        idx = idx + 360;
    end
    set(h, 'Color', cmap(idx,:)) % before, change each quiver command to h = quiver(...)
    hold on
    
    for z=1:length(iso)
        h = quiver3(iso(z).x,iso(z).y,0,iso(z).Mx(i)*norm,iso(z).My(i)*norm,iso(z).Mz(i)*norm,'linewidth',2,'LineStyle','-');
        alpha = atan2(iso(z).My(i)*norm, iso(z).Mx(i)*norm);
        idx = ceil(rad2deg(alpha));
        if idx < 1 % negative indices and 0 are invalid
            idx = idx + 360;
        end
        set(h, 'Color', cmap(idx,:)) % before, change each quiver command to h = quiver(...)
    end
    
    grid on
    box on
    hold off
    xlim([min(x)-norm*norm_scale max(x)+norm*norm_scale]);
    %ylim([-norm*norm_scale norm*norm_scale]);
    ylim([min(x)-norm*norm_scale max(x)+norm*norm_scale]);
    zlim([min(x)-norm*norm_scale max(x)+norm*norm_scale]);
    xlabel("$M_{x'}$", "interpreter", "latex", "fontsize", 10)
    ylabel("$M_{y'}$", "interpreter", "latex", "fontsize", 10)
    zlabel("$M_{z'}$", "interpreter", "latex", "fontsize", 10)
    
    h1 = subplot(2,2,3:4);
    title(['Time: ', num2str(Time(i)*1000), ' ms'])
    %     plot(h1,Time(1:i),Mx(1:i),'linewidth',1,'LineStyle','-','Color','r')
    hold on
    plot(h1,Time(1:i)*1000,My(1:i),'linewidth',1,'LineStyle','-','Color','b')
    %     plot(h1,Time(1:i),Mz(1:i),'linewidth',1,'LineStyle','-','Color','k')
    grid on
    box on
    hold off
    if Time(i)<(0+.32+5.12)/1000
        xlim(h1,[0 (0+.32+5.12)]);
    else
        xlim(h1,[0 max(Time)*1000]);
    end
    ylim(h1,[-0.2 1.2]);
    xlabel(h1,"Time (ms)", "interpreter", "latex", "fontsize", 10)
    ylabel(h1,"Signal", "interpreter", "latex", "fontsize", 10)
    
    h1 = subplot(2,2,2);
    plot(h1,iso(1).My(1:i),iso(1).Mz(1:i),'linewidth',1,'LineStyle','--','Color','r')
    hold on
    quiver(h1,0,0,iso(1).My(i),iso(1).Mz(i),'linewidth',2,'LineStyle','-','Color','r')
    grid on
    box on
    hold off
    legend(h1,'Path','Now', "interpreter", "latex", "fontsize", 10)
    xlim(h1,[-0.2 1.2]);
    ylim(h1,[-0.2 1.2]);
    xlabel(h1,"$M_{y'}$", "interpreter", "latex", "fontsize", 10)
    ylabel(h1,"$M_{z'}$", "interpreter", "latex", "fontsize", 10)
    pause(0.01)
end





