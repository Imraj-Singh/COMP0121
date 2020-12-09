%% Problem 1: Experiment 2 gradient relaxation
% Author: Imraj Singh 08/11/2020

clc
clear
close all

% Add the path that holds functions
addpath 'C:\Users\Imraj Singh\Documents\UCL\Comp_MRI\COMP0121\Coursework_2\Functions'

%% Define coordinate position of each isochromat

% Length of segment
L = 20/1000;

% Sampling time
Ts = 5.12;

% Sampling time
Gradient = -4.6;

% Number of isochromats x
Numx = 101;

% Number of isochromats y
Numy = 101;

% Spatial coordinates of isochromats (m)
x1 = linspace(-L/2,0,Numx);
y1 = linspace(0,L/2,Numy);

% Create mesh grid of x and y coordinates
[X1, Y1] = meshgrid(x1,y1);

% Spatial coordinates of isochromats (m)
x2 = linspace(0,L/2,Numx);
y2 = linspace(-L/2,0,Numy);

% Create mesh grid of x and y coordinates
[X2, Y2] = meshgrid(x2,y2);

% Flatten the array into a vector to make manipulation easier
X = [X1(:);X2(:)];
Y = [Y1(:);Y2(:)];

%% Set NMR parameters and initialise isochromat structures

% Gyromagnetic ratio of hyrodrogen proton
gamma = 2.68*10^8;

% Initalise the isochromats, setting the X, Y cooridates and instantiating
% the Mx, My, Mz fields of the structures
iso = initialise_iso(X, Y);

% B0
B0 = 1.5;

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
%blocks.N = [1 65 129 257];
blocks.N = [1 11 101 201];

% Time over which the block occurs
blocks.t = [0 .32 Ts/2 Ts]/1000;

% The gradient of the field dB/dx for each block
blocks.Gx = [0 0 Gradient -Gradient]/1000;

% The gradient of the field dB/dz for each block
blocks.Gy = [0 0 Gradient -Gradient]/1000;

% The relaxation time only when spins are relaxation
blocks.relax = [0 0 Ts/2 Ts]/1000;

% Initialise the Mx, My, Mz values that stores the sum of all the spins
% magnetisation
Mx = zeros(sum(blocks.N),1);
My = zeros(sum(blocks.N),1);
Mz = zeros(sum(blocks.N),1);

map = ones(Numx*2-1,Numy*2-1,sum(blocks.N))*361;
map(1,1,:) = 0;

% Calculate magnetisation for all the spins
for i=1:length(iso)
    % Initialise spin
    iso(i) = initial_dir(iso(i),'z',blocks.N(1));
    
    % RF pulse
    iso(i) = rf_pulse(iso(i), 'x', pi/2, blocks.N(2), blocks.t(2));
    
    % Free relaxation with gradient
    iso(i) = free_relaxation(iso(i), blocks.N(3), blocks.t(3), ...
        blocks.Gx(3), blocks.Gy(3), blocks.relax(1:3), gamma, T1, T2, R_FoR);
    
    % Free relaxation with negative gradient
    iso(i) = free_relaxation(iso(i), blocks.N(4), blocks.t(4), ...
        blocks.Gx(4), blocks.Gy(4), blocks.relax(1:4), gamma, T1, T2, R_FoR);
    
    ix = round((iso(i).x+L/2)/(x1(2)-x1(1)) + 1);
    iy = round((iso(i).y+L/2)/(y1(2)-y1(1)) + 1);
    
    for k=1:sum(blocks.N)
        map(ix,iy,k) = ceil(mod(atan2(iso(i).My(k), iso(i).Mx(k)),2*pi)/(pi)*180);
    end
    
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

colormap([1 1 1; hsv(360); 0 0 0]);
for i=1:length(Time)
imagesc(linspace(-L/2,L/2,Numx*2-1),linspace(-L/2,L/2,Numy*2-1),map(:,:,i))
set(gca,'YDir','normal')
pause(0.1)

end


%% Animation module

% Normalisation constants used for plotting the quiver3s
norm = 0.001;
norm_scale = 2;

% Create circular colour map
cmap = colormap(hsv(360));

% Set the video writing name and location
video = VideoWriter(['P1_E3', '.mp4'], 'MPEG-4');

% Set the frame rate of video
frameRate = 25;
video.set('FrameRate', frameRate);

% Open the video
video.open();

% Open a figure that has frame recorded
h = figure;

for i=1:length(Time)
    
    if Time(i)>=0 && Time(i)<=sum(blocks.t(1:1))
    sgtitle(['Flipping the spins, time: ', num2str(Time(i)*1000), ' ms'], "interpreter", "latex", "fontsize", 15)
    elseif Time(i)>=sum(blocks.t(1:1)) && Time(i)<=sum(blocks.t(1:2))
    sgtitle(['Flipping the spins, time: ', num2str(Time(i)*1000), ' ms'], "interpreter", "latex", "fontsize", 15)
    elseif Time(i)>sum(blocks.t(1:2)) && Time(i)<=sum(blocks.t(1:3))
    sgtitle(['Dephasing lobe, time (after flip): ', num2str((Time(i)-blocks.t(2))*1000), ' ms'], "interpreter", "latex", "fontsize", 15)
    elseif Time(i)>sum(blocks.t(1:3)) && Time(i)<=sum(blocks.t(1:4))
    sgtitle(['Recording signal, time (after dephase): ', num2str((Time(i)-blocks.t(2))*1000), ' ms'], "interpreter", "latex", "fontsize", 15)
    end
    
    % First subplot - quiver of spin vectors in space
    subplot(2,2,1)
    % First element of vector of iso structs
    h1 = quiver3(iso(1).x,iso(1).y,0,iso(1).Mx(i)*norm,iso(1).My(i)*norm,iso(1).Mz(i)*norm,'linewidth',2,'LineStyle','-','Color','k');
    % Calculate the angle of the spin from the positive y'-axis clockwise
    alpha = atan2(iso(1).My(i), iso(1).Mx(i));
    % Convert the angle to degrees and round
    idx = ceil(rad2deg(alpha));
    % If the degree are negative/zero then invalid and add 360 degrees
    if idx < 1
        idx = idx + 360;
    end
    % Set the colour of the quiver according to the colour map
    set(h1, 'Color', cmap(idx,:))
    hold on
    
    for z=2:length(iso)
        % z'th element of vector of iso structs
        h1 = quiver3(iso(z).x,iso(z).y,0,iso(z).Mx(i)*norm,iso(z).My(i)*norm,iso(z).Mz(i)*norm,'linewidth',2,'LineStyle','-','Color','k');
        % Calculate the angle of the spin from the positive y'-axis clockwise
        alpha = atan2(iso(z).My(i), iso(z).Mx(i));
        % Convert the angle to degrees and round
        idx = ceil(rad2deg(alpha));
        % If the degree are negative/zero then invalid and add 360 degrees
        if idx < 1
            idx = idx + 360;
        end
        % Set the colour of the quiver according to the colour map
        set(h1, 'Color', cmap(idx,:))
        hold on
    end
    
    grid on
    box on
    hold off
    % Set the axes limits
    xlim([min(x)-norm*norm_scale max(x)+norm*norm_scale]);
    ylim([min(x)-norm*norm_scale max(x)+norm*norm_scale]);
    zlim([min(x)-norm*norm_scale max(x)+norm*norm_scale]);
    % Change to 2D view after flip
    if Time(i) > blocks.t(2)
        view(2)
    end
    % Set the axes labels
    xlabel("$x'$", "interpreter", "latex", "fontsize", 10)
    ylabel("$y'$", "interpreter", "latex", "fontsize", 10)
    zlabel("$z'$", "interpreter", "latex", "fontsize", 10)
    
    % Second subplot - Magnetisation vs Time
    subplot(2,2,2);
    % Plot the full magnetisation
    plot(Time(1:i)*1000,Mx(1:i),'linewidth',1,'LineStyle','-','Color','r')
    hold on
    plot(Time(1:i)*1000,My(1:i),'linewidth',1,'LineStyle','-','Color','b')
    % plot the signal points
    if Time(i)<sum(blocks.t(1:4)) &&  Time(i)>sum(blocks.t(1:3))
        plot(Time(sum(blocks.N(1:3)):i)*1000,Mx(sum(blocks.N(1:3)):i),'.','MarkerSize',8,'Color','k')
        plot(Time(sum(blocks.N(1:3)):i)*1000,My(sum(blocks.N(1:3)):i),'.','MarkerSize',8,'Color','k')
    end
    grid on
    box on
    hold off
    % Set the legend
    legend("$M_{x'}$","$M_{y'}$", "interpreter", "latex", "fontsize", 10,'Location','southeast')
    % Set the axes limits
    xlim([0 max(Time)*1.2*1000]);
    ylim([-1 1]*1.2);
    % Set the axes labels
    xlabel("Time (ms)", "interpreter", "latex", "fontsize", 10)
    ylabel("Magnetisation", "interpreter", "latex", "fontsize", 10)
    
    % Third subplot - K space values
    subplot(2,2,3);
    plot(kx(1:i)/1000,ky(1:i)/1000,'linewidth',1,'LineStyle','-','Color','r')
    hold on
    % Current K space value
    plot(kx(i)/1000,ky(i)/1000,'.','MarkerSize',8,'Color','k')
    grid on
    box on
    hold off
    % Set the axes limits
    xlim([-1.2 1.2]);
    ylim([-1.2 1.2]);
    % Set the axes labels
    xlabel("$k_x$ (1/mm)", "interpreter", "latex", "fontsize", 10)
    ylabel("$k_y$ (1/mm)", "interpreter", "latex", "fontsize", 10)
    
    % Fourth subplot - Time vs K space
    subplot(2,2,4);
    % X component
    plot(Time(1:i)*1000,kx(1:i)/1000,'linewidth',1,'LineStyle','-','Color','r')
    hold on
    % Y component
    plot(Time(1:i)*1000,ky(1:i)/1000,'linewidth',1,'LineStyle','-','Color','b')
    grid on
    box on
    hold off
    % Set the legend
    legend("$k_x$","$k_y$", "interpreter", "latex", "fontsize", 10,'Location','northwest')
    % Set the axes limits
    xlim([0 max(Time)*1.2*1000]);
    ylim([-1.2 1.2])
    % Set the axes labels
    xlabel("Time (ms)", "interpreter", "latex", "fontsize", 10)
    ylabel("Spatial frequency (1/mm)", "interpreter", "latex", "fontsize", 10)
    
    % Set frame to add to video
    frame = getframe(h);
    video.writeVideo(frame);
end

% Close the video
video.close();
