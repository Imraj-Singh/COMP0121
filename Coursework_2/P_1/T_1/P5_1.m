%% Free procession with Dephasing
% Author: Imraj Singh 03/11/2020

clc
clear
close all

addpath(genpath('.../Functions'))
% Number of isochromats x
Numx = 11;
% Number of isochromats y
Numy = 11;
% Spatial coordinates of isochromats (m)
x = linspace(-10,10,Numx)/1000;
y = linspace(-10,10,Numy)/1000;

[X, Y] = meshgrid(x,y);

X = X(:);
Y = Y(:);


% Gyromagnetic ratio
gamma = 2.68*10^8;
% Initalise the isochromats
iso = initalise_iso(X, Y);
% B0
B0 = 3;
% Initalise the B0
iso = initalise_B0(iso, B0);
% Larmor frequency of rotating frame of reference
R_FoR = B0*gamma;
%
T2 = 1;
T1 = 2;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orientate -> rf pulse -> free precession
% [N, t, Gx, Gy, relax]
blocks.N = [3 11 11 11 11 11];
blocks.t = [0 3.2 5.12 5.12 5.12 5.12]/1000;
blocks.Gx = [0 0 4.6 -4.6 4.6 -4.6]/1000;
blocks.Gy = [0 0 0 0 0 0]/1000;
blocks.relax = [0 0 5.12 5.12 5.12 5.12]/1000;

Mx = zeros(sum(blocks.N),1);
My = zeros(sum(blocks.N),1);
Mz = zeros(sum(blocks.N),1);

for i=1:length(iso)
        iso(i) = initial_dir(iso(i),'z',blocks.N(1));
        
        iso(i) = rf_pulse(iso(i), 'x', pi/2, blocks.N(2), blocks.t(2));
        
        iso(i) = free_relaxation(iso(i), blocks.N(3), blocks.t(3), blocks.Gx(3), blocks.Gy(3), blocks.relax, gamma, T1, T2, R_FoR);
        
        iso(i) = free_relaxation(iso(i), blocks.N(4), blocks.t(4), blocks.Gx(4), blocks.Gy(4), blocks.relax, gamma, T1, T2, R_FoR);
        
        iso(i) = free_relaxation(iso(i), blocks.N(5), blocks.t(5), blocks.Gx(5), blocks.Gy(5), blocks.relax, gamma, T1, T2, R_FoR);
        
        iso(i) = free_relaxation(iso(i), blocks.N(6), blocks.t(6), blocks.Gx(6), blocks.Gy(6), blocks.relax, gamma, T1, T2, R_FoR);
%         for k=1:sum(blocks.N)
%             map(i,k) = atan2(iso(i).My(k), iso(i).Mx(k));
%         end
        Mx = Mx + iso(i).Mx;
        My = My + iso(i).My;
        Mz = Mz + iso(i).Mz;
end

Mx = Mx/(Numx*Numy);
My = My/(Numx*Numy);
Mz = Mz/(Numx*Numy);

Time = linspace(0, blocks.t(1), blocks.N(1));
for i = 2:length(blocks.t)
   Time = [Time, linspace(sum(blocks.t(1:(i-1))),sum(blocks.t(1:i)),blocks.N(i))];
end


%% Animation module
norm = 0.001;
norm_scale = 2;

Gx_t = [];
Gy_t = [];
kx = [];
ky = [];

[Gx_t, Gy_t, kx, ky] = unroll_plotting(blocks.N,blocks.Gx,blocks.Gy,gamma,Time);

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
    
    view(2)
    grid on
    box on
    hold off
    xlim([min(x)-norm*norm_scale max(x)+norm*norm_scale]);
    %ylim([-norm*norm_scale norm*norm_scale]);
    ylim([min(x)-norm*norm_scale max(x)+norm*norm_scale]);
    zlim([0 norm]);
    xlabel("$M_{x'}$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_{y'}$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_{z'}$", "interpreter", "latex", "fontsize", 15)
    
    h1 = subplot(2,2,2);
    plot(h1,Time(1:i),Mx(1:i),'linewidth',1,'LineStyle','-','Color','r')
    hold on
    plot(h1,Time(1:i),My(1:i),'linewidth',1,'LineStyle','-','Color','b')
    grid on
    box on
    hold off
    legend(h1,'x','y')
    xlim(h1,[0 max(Time)*1.2]);
    ylim(h1,[-1 1]*1.2);
    xlabel(h1,"Time (s)", "interpreter", "latex", "fontsize", 15)
    ylabel(h1,"Signal", "interpreter", "latex", "fontsize", 15)
    
    h1 = subplot(2,2,3);
    plot(h1,kx(1:i),ky(1:i),'linewidth',1,'LineStyle','-','Color','r')
    hold on
%     plot(h1,Time(1:i),ky(1:i),'linewidth',1,'LineStyle','-','Color','b')
    grid on
    box on
    hold off
    xlim(h1,[min([min(kx) min(ky)])*1.2 max([max(kx) max(ky)])*1.2]);
    ylim(h1,[min([min(kx) min(ky)])*1.2 max([max(kx) max(ky)])*1.2]);
    xlabel(h1,"Time (s)", "interpreter", "latex", "fontsize", 15)
    ylabel(h1,"Signal", "interpreter", "latex", "fontsize", 15)
    
    h1 = subplot(2,2,4);
    plot(h1,Time(1:i),kx(1:i),'linewidth',1,'LineStyle','-','Color','r')
    hold on
    plot(h1,Time(1:i),ky(1:i),'linewidth',1,'LineStyle','-','Color','b')
    grid on
    box on
    hold off
    legend(h1,'x','y')
    xlim(h1,[0 max(Time)*1.2]);
    ylim(h1, [min([min(kx) min(ky)])*1.2 max([max(kx) max(ky)])*1.2])
    xlabel(h1,"Time (s)", "interpreter", "latex", "fontsize", 15)
    ylabel(h1,"Signal", "interpreter", "latex", "fontsize", 15)
    pause(0.01)
end





