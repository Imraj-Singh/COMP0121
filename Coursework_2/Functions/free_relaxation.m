function iso = free_relaxation(iso,N,T,Gx,Gy,relax, gamma, T1, T2,R_FoR)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

omega = gamma*iso.B0 + (Gx*iso.x + Gy*iso.y) * gamma / (2*pi);
omega = omega - R_FoR;
t_Relax = linspace(sum(relax(1:(end-1))),sum(relax),N)';

M_z_begin = (iso.Mz(end)-(1-exp(-t_Relax(1)/T1)))/exp(-t_Relax(1)/T1);
t = linspace(0,T,N)';


M_x = (iso.Mx(end).*cos(omega*t) + iso.My(end).*sin(omega*t)).*exp(-t_Relax/T2);
M_y = (iso.My(end).*cos(omega*t) - iso.Mx(end).*sin(omega*t)).*exp(-t_Relax/T2);
M_z = M_z_begin*exp(-t_Relax/T1) + (1-exp(-t_Relax./T1));

iso.Mx = [iso.Mx; M_x];
iso.My = [iso.My; M_y];
iso.Mz = [iso.Mz; M_z];
end

