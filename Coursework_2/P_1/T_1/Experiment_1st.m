%% Problem 1: Experiment 1 free precession
% Author: Imraj Singh 08/11/2020

clc
clear
close all

% Add the path that holds functions
addpath(genpath('..\Functions'))

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
iso = initialise_iso(X, Y);

% B0
B0 = 3;

% Initalise the B0 this is done such that arbitrary perturbation can be
% added to the B0 field to simulate field inhomogenities
iso = initialise_B0(iso, B0);

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
blocks.N = [11 101 257 151];

% Time over which the block occurs
blocks.t = [0 .32 5.12 4994.88]/1000;

% The gradient of the field dB/dx for each block
blocks.Gx = [0 0 0 0]/1000;

% The gradient of the field dB/dz for each block
blocks.Gy = [0 0 0 0]/1000;

% The relaxation time only when spins are relaxation
blocks.relax = [0 0 5.12 4994.88]/1000;

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

% Normalisation constants used for plotting the quiver3s
norm = 0.001;
norm_scale = 2;

% Create circular colour map
cmap = colormap(hsv(360));

% Set the video writing name and location
video = VideoWriter(['P1_E1', '.mp4'], 'MPEG-4');

% Set the frame rate of video
frameRate = 50;
video.set('FrameRate', frameRate);

% Open the video
video.open();

% Open a figure that has frame recorded
h = figure;

for i=1:length(Time)
    % First subplot
    subplot(2,2,1)
    % First element of vector of iso structs
    quiv = quiver3(iso(1).x,iso(1).y,0,iso(1).Mx(i)*norm,iso(1).My(i)*norm,iso(1).Mz(i)*norm,'linewidth',2,'LineStyle','-','Color','k');
    % Calculate the angle of the spin from the positive y'-axis clockwise
    alpha = atan2(iso(1).My(i), iso(1).Mx(i));
    % Convert the angle to degrees and round
    idx = ceil(rad2deg(alpha));
    % If the degree are negative/zero then invalid and add 360 degrees
    if idx < 1
        idx = idx + 360;
    end
    % Set the colour of the quiver according to the colour map
    set(quiv, 'Color', cmap(idx,:))
    hold on
    
    for z=2:length(iso)
        % z'th element of vector of iso structs
        quiv = quiver3(iso(z).x,iso(z).y,0,iso(z).Mx(i)*norm,iso(z).My(i)*norm,iso(z).Mz(i)*norm,'linewidth',2,'LineStyle','-','Color','k');
        % Calculate the angle of the spin from the positive y'-axis clockwise
        alpha = atan2(iso(1).My(i), iso(1).Mx(i));
        % Convert the angle to degrees and round
        idx = ceil(rad2deg(alpha));
        % If the degree are negative/zero then invalid and add 360 degrees
        if idx < 1
            idx = idx + 360;
        end
        % Set the colour of the quiver according to the colour map
        set(quiv, 'Color', cmap(idx,:))
        hold on
    end
    
    grid on
    box on
    hold off
    % Set the axes limits
    xlim([min(x)-norm*norm_scale max(x)+norm*norm_scale]);
    ylim([min(x)-norm*norm_scale max(x)+norm*norm_scale]);
    zlim([min(x)-norm*norm_scale max(x)+norm*norm_scale]);
    % Set the axes labels
    xlabel("$x'$", "interpreter", "latex", "fontsize", 15)
    ylabel("$y'$", "interpreter", "latex", "fontsize", 15)
    zlabel("$z'$", "interpreter", "latex", "fontsize", 15)
    
    % Second subplot
    subplot(2,2,2);
    % Plot the path of the magnetisation vector
    plot(iso(1).My(1:i),iso(1).Mz(1:i),'linewidth',1,'LineStyle','--','Color','r')
    hold on
    % Plot the magnetisation vector
    quiver(0,0,iso(1).My(i),iso(1).Mz(i),'linewidth',2,'LineStyle','-','Color','r')
    grid on
    box on
    hold off
    % Set the axes limits
    xlim([-0.2 1.2]);
    ylim([-0.2 1.2]);
    % Set the axes labels
    xlabel("$M_{y'}$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_{z'}$", "interpreter", "latex", "fontsize", 15)
    
    % Last subplot, double length
    subplot(2,2,3:4);
    % Set the title based on what part of the pulse sequence is animated
    if Time(i)<=(blocks.t(2))
    title(['Flipping the spins, time: ', num2str(Time(i)*1000), ' ms'], "interpreter", "latex", "fontsize", 15)
    elseif Time(i)<=(blocks.t(2)+blocks.t(3))
    title(['Recording signal, time (after flip): ', num2str((Time(i)-blocks.t(2))*1000), ' ms'], "interpreter", "latex", "fontsize", 15)
    elseif Time(i)>(blocks.t(2)+blocks.t(3))
    title(['Full relaxation, time (after flip): ', num2str(ceil((Time(i)-blocks.t(2))*1000)), ' ms'], "interpreter", "latex", "fontsize", 15)
    end
    hold on
    % Plot the y' and z' magnetisation as these are the most important 
    plot(Time(1:i)*1000,My(1:i),'linewidth',2,'LineStyle','-','Color','b')
    plot(Time(1:i)*1000,Mz(1:i),'linewidth',2,'LineStyle','-','Color','r')
    grid on
    box on
    hold off
    % Set two different scales to investigate the full relaxation
    if Time(i)<(blocks.t(2)+blocks.t(3))
        xlim([0 (blocks.t(2)+blocks.t(3))*1000]);
    else
        xlim([0 max(Time)*1000]);
    end
    % Add legends, limits, labels
    legend("$M_{y'}$","$M_{z'}$", "interpreter", "latex", "fontsize", 15)
    ylim([-0.2 1.2]);
    xlabel("Time (ms)", "interpreter", "latex", "fontsize", 15)
    ylabel("Magnetisation", "interpreter", "latex", "fontsize", 15)
    
    % Set frame to add to video
    frame = getframe(h);
    video.writeVideo(frame);
end

% Close the video
video.close();



