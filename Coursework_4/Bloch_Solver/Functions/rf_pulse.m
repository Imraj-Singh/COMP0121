function iso = rf_pulse(iso, Axis, Angle, N, T)
%FREE_RELAXATION by Imraj Singh
% Inputs
%   iso, Axis, Angle, N, T
% Outputs
%   iso

% Error handling
if any((Axis~='x')&(Axis~='y')&(Axis~='z'))
    error('Needs be an axis of either "x","y" or "z".')
end
if (~isa(Angle,'double'))
    error('Angle needs to be a double.')
end

% Create time vector
t = linspace(0,T,N)';

% Calculate
if (Axis=='x')
    M_x = iso.Mx(end)*ones(N,1);
    M_y = iso.My(end).*cos(Angle/T*t) + iso.Mz(end).*sin(Angle/T*t);
    M_z = iso.Mz(end).*cos(Angle/T*t) - iso.My(end).*sin(Angle/T*t);
elseif (Axis=='y')
    M_x = iso.Mx(end).*cos(Angle/T*t) - iso.Mz(end).*sin(Angle/T*t);
    M_y = iso.My(end).*ones(N,1);
    M_z = iso.Mz(end).*cos(Angle/T*t) + iso.Mx(end).*sin(Angle/T*t);
elseif (Axis=='z')
    M_x = iso.Mx(end).*cos(Angle/T*t) + iso.My(end).*sin(Angle/T*t);
    M_y = iso.My(end).*cos(Angle/T*t) - iso.Mx(end).*sin(Angle/T*t);
    M_z = iso.Mz(end)*ones(N,1);
end

% Add it appropriate fields
iso.Mx = [iso.Mx; M_x];
iso.My = [iso.My; M_y];
iso.Mz = [iso.Mz; M_z];
end

