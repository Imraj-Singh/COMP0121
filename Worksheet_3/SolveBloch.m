function SolveBloch(time)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
timeV = linspace(0, time, 101);

Mz = zeros([1,length(timeV)]);
Mx = zeros([1,length(timeV)]);
My = zeros([1,length(timeV)]);

Mz(1) = 0;
Mx(1) = 1;
My(1) = 1;

gamma = 2.68*10^8;
b0 = 7;
omega0 = gamma * b0;
M0 = 100;
T1 = 1;
T2 = T1/5;

for i = 1:(length(timeV)-1)
    Mz(i+1) = Mz(i) + (M0 - Mz(i))/T1;
    Mx(i+1) = omega0 * My(i) + Mx(i)/T2;
    My(i+1) = omega0 * Mx(i) + My(i)/T2;
end
figure
plot3(Mx,My,Mz)

end