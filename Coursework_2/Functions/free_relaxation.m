function iso = free_relaxation(iso,N,T,Gx,Gy,relax, gamma, T1, T2,R_FoR)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

omega = (iso.B0 + (Gx*iso.x + Gy*iso.y)) * gamma;
omega = omega - R_FoR;
t_Relax = linspace(sum(relax(1:(end-1))),sum(relax),N)';

M_x_begin = iso.Mx(end)/exp(-t_Relax(1)/T2);
M_y_begin = iso.My(end)/exp(-t_Relax(1)/T2);
M_z_begin = (iso.Mz(end)-(1-exp(-t_Relax(1)/T1)))/exp(-t_Relax(1)/T1);
t = linspace(0,T,N)';


M_x = (M_x_begin.*cos(omega*t) + M_y_begin.*sin(omega*t)).*exp(-t_Relax/T2);
M_y = (M_y_begin.*cos(omega*t) - M_x_begin.*sin(omega*t)).*exp(-t_Relax/T2);
M_z = M_z_begin*exp(-t_Relax/T1) + (1-exp(-t_Relax./T1));

iso.Mx = [iso.Mx; M_x];
iso.My = [iso.My; M_y];
iso.Mz = [iso.Mz; M_z];
end

