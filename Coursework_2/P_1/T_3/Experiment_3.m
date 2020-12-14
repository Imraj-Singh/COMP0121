%% Problem 1: Experiment 2 gradient relaxation
% Author: Imraj Singh 08/11/2020

clc
clear
close all

% Add the path that holds functions
addpath 'C:\Users\Imraj Singh\Documents\UCL\Comp_MRI\COMP0121\Coursework_2\Functions'

%% Define coordinate position of each isochromat

% Number of isochromats x
Numx = 100;

% Number of isochromats y
Numy = 101;

% Length of segment
L = 20/1000;

% Spatial spacing
dx = L/Numx;

% Spatial coordinates of isochromats (m)
x = linspace(-L/2 + dx/2,L/2 - dx/2,Numx);
y = linspace(-L/2,L/2,Numy);

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

% Time of sample
Ts = 5.12;

% Each experiment begins with initialising the direction of the spins then
% by flipping the spins with a pi/2 hard RF pulse along the x' axis over
% 0.32ms

% [Set the number of points to be calculated minimum of 2 for each block
blocks.N = [1 65 129 257];

% Time over which the block occurs
blocks.t = [0 .32 Ts/2 Ts]/1000;

% The gradient of the field dB/dx for each block
blocks.Gy = [0 0 -4.6 4.6]/1000;

% The gradient of the field dB/dz for each block
blocks.Gx = [0 0 0 0]/1000;

% The relaxation time only when spins are relaxation
blocks.relax = [0 0 Ts/2 Ts]/1000;

% Reset k-space
blocks.resetk = [1 0 0 0];

% Where the signal is collected
blocks.signal = [0 0 0 1];

% Initialise the Mx, My, Mz values that stores the sum of all the spins
% magnetisation
Mx = zeros(sum(blocks.N),1);
My = zeros(sum(blocks.N),1);
Mz = zeros(sum(blocks.N),1);

% Subset definition
Sub_num = 20;
Sub_set = linspace(-L/2 + (L/Sub_num)/2,L/2 - (L/Sub_num)/2,Sub_num);

% Colour map ensuring that only isochromats have colour
map = ones(Numx,Numx,sum(blocks.N))*-1;

z = 1;
% Calculate magnetisation for all the spins
for i=1:length(iso)
    % Initialise spin
    iso(i) = initial_dir(iso(i),'z',blocks.N(1));
    
    % RF pulse
    iso(i) = rf_pulse(iso(i), 'x', pi/2, blocks.N(2), blocks.t(2));
    
    % Free relaxation with gradient
    iso(i) = free_relaxation(iso(i), blocks.N(3), blocks.t(3), ...
        blocks.Gx(3), blocks.Gy(3), blocks.relax(1:3), gamma, T1, T2, R_FoR);
    
    
    % Free relaxation with gradient
    iso(i) = free_relaxation(iso(i), blocks.N(4), blocks.t(4), ...
        blocks.Gx(4), blocks.Gy(4), blocks.relax(1:4), gamma, T1, T2, R_FoR);
    
    
    % Calculate the index's for the map array dependent on spatial location
    ix = round((iso(i).x+L/2)/(x(2)-x(1)) + 1);
    iy = round((iso(i).y+L/2)/(y(2)-y(1)) + 1);
    
    % Calculate the angles on the x-y plane the spin vector for each
    % isochromat corresponding to a pixel, additionally we round up with
    % each pixel as there are only 360 unique colours in the colour map,
    % also we made sure that the calculated angles fall between 1->360
    for k=1:sum(blocks.N)
        map(iy,ix,k) = ceil(mod(atan2(iso(i).My(k), iso(i).Mx(k)),2*pi)/(pi)*180);
    end
    
    
    % Find the subset of isos that we specify that can select certain
    % isochromats within the whole domain which we calculated for the
    % imagesc for better visualisation
    if  ismembertol(iso(i).x,Sub_set) && ismembertol(iso(i).y,0)
        % Add the magnetisations to the total  magnetisation
        Mx = Mx + iso(i).Mx;
        My = My + iso(i).My;
        Mz = Mz + iso(i).Mz;
        Sub_set_iso(z) = iso(i);
        z=z+1;
    end
    
    % Clear the arrays as they take up a lot of memory
    iso(i).Mx = [];
    iso(i).My = [];
    iso(i).Mz = [];
end

% Magnetisation now normalised across all the spins
Mx = Mx/length(Sub_set_iso);
My = My/length(Sub_set_iso);
Mz = Mz/length(Sub_set_iso);

% Calculate the time vector
Time = unroll_time(blocks.N,blocks.t);

% Calculate the Gx, Gy, kx, ky vectors in time
[Gx_t, Gy_t, kx, ky] = unroll_plotting(blocks.N,blocks.Gx,blocks.Gy,blocks.resetk,gamma,Time);

% Binary filter used to just plot points where the signal is taken
s_filter = unroll_signal(blocks.N, blocks.signal);

%% Animation module

% Normalisation constants used for plotting the quiver3s
norm = 1;
norm_scale = 2;

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
    colormap([1 1 1; hsv(360)]);
    
    if Time(i)>=0 && Time(i)<sum(blocks.t(1:2))
        imagesc(linspace(-L/2,L/2,Numx*2)*1000,linspace(-L/2,L/2,Numx*2-1)*1000,ones(length(map(:,:,1)))*-1,[-1 360])
    elseif Time(i)>=sum(blocks.t(1:2))
        imagesc(linspace(-L/2,L/2,Numx*2)*1000,linspace(-L/2,L/2,Numx*2-1)*1000,map(:,:,i),[-1 360])
    end
    hold on
    %First element of vector of iso structs
    quiver3(Sub_set_iso(1).x*1000,Sub_set_iso(1).y*1000,0,Sub_set_iso(1).Mx(i)*norm,Sub_set_iso(1).My(i)*norm,Sub_set_iso(1).Mz(i)*norm,'linewidth',2,'LineStyle','-','Color','k');
    scatter(Sub_set_iso(1).x*1000,Sub_set_iso(1).y*1000, 10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'LineWidth',1);
    for z=2:length(Sub_set_iso)
        % z'th element of vector of iso structs
        quiver3(Sub_set_iso(z).x*1000,Sub_set_iso(z).y*1000,0,Sub_set_iso(z).Mx(i)*norm,Sub_set_iso(z).My(i)*norm,Sub_set_iso(z).Mz(i)*norm,'linewidth',2,'LineStyle','-','Color','k');
        scatter(Sub_set_iso(z).x*1000,Sub_set_iso(z).y*1000, 10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'LineWidth',1);
    end
    set(gca,'YDir','normal');
    grid on
    box on
    hold off
    % Set the axes limits
    xlim([-L/2 L/2]*1000);
    ylim([-L/2 L/2]*1000);
    % Change to 2D view after flip
    if Time(i) > blocks.t(2)
        view(2)
    end
    % Set the axes labels
    xlabel("$x$ (mm)", "interpreter", "latex", "fontsize", 10)
    ylabel("$y$ (mm)", "interpreter", "latex", "fontsize", 10)
    
    % Second subplot - Magnetisation vs Time
    subplot(2,2,2);
    % Plot the full magnetisation
    plot(Time(1:i)*1000,Mx(1:i),'linewidth',1,'LineStyle','-','Color','r')
    hold on
    plot(Time(1:i)*1000,My(1:i),'linewidth',1,'LineStyle','-','Color','b')
    % plot the signal points
    plot(Time(find(s_filter(1:i)))*1000,My(find(s_filter(1:i))),'.','MarkerSize',5,'Color','k')
    grid on
    box on
    hold off
    % Set the legend
    legend("$M_{x'}$","$M_{y'}$", "interpreter", "latex", "fontsize", 10,'Location','southeast')
    % Set the axes limits
    xlim([0 max(Time)*1000]);
    ylim([-1 1]*1.2);
    % Set the axes labels
    xlabel("Time (ms)", "interpreter", "latex", "fontsize", 10)
    ylabel("Magnetisation", "interpreter", "latex", "fontsize", 10)
    
    % Third subplot - K space values
    subplot(2,2,3);
    % plot the current location
    plot(kx(i)/1000,ky(i)/1000,'x','MarkerSize',5,'Color','k')
    hold on
    % plot the signal points
    if Time(i)>=sum(blocks.t(1:2))
        plot(kx(find(s_filter(1:i)))/1000,ky(find(s_filter(1:i)))/1000,'.','MarkerSize',8,'Color','k')
    end
    grid on
    box on
    hold off
    % Set the axes limits
    xlim([-1.2 1.2]);
    ylim([-1.2 .2]);
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
    xlim([0 max(Time)*1000]);
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
