function SE = calcSE(T,B, Gamma)
%Function to calculate the approximate spin excess
%   Inputs
%   T = Temperature of the medium
%   B = Magnetic field strength
%   Gamma = Gyromagnetic ratio of medium
%
%   Output
%   SE = Spin excess

SE = T + B + Gamma;
end

