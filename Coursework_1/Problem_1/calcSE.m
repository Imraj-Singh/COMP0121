function SE = calcSE(T,B, Gamma)
%Function to calculate the approximate spin excess
%   Inputs
%   T = Temperature of the medium
%   B = Magnetic field strength
%   Gamma = Gyromagnetic ratio of medium
%
%   Output
%   SE = Spin excess

% Parameters

% Number of spins calculated from typical voxel volume 2*2*5mm = 0.02
% Avogadro's number = 6.02*10^23
% One voxel has 2 x 6.02 x10^23 x 0.02 / 18 protons
N = 1.338 * 10^21;

% Reduced planks constant
hbar = (6.62607015*10^(-34))/(2*pi);

% Boltzman constant
k = 1.38064852 * 10^(-23);

SE = N*(hbar*B*Gamma)/(2*k*T);
end

