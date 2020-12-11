function [Gx_t, Gy_t, kx, ky] = unroll_plotting(N,Gx,Gy,resetk,gamma,Time)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Gx_t = [];
Gy_t = [];
kx = zeros(1,sum(N));
ky = zeros(1,sum(N));


for i=1:length(N)
    Gx_t = [Gx_t, ones(1,N(i))*Gx(i)];
    Gy_t = [Gy_t, ones(1,N(i))*Gy(i)];
end

for i=2:length(N)
    if resetk(i)==1
        kx(sum(N(1:i-1)):sum(N(1:i))) = linspace(0,0,N(i));
        ky(sum(N(1:i-1)):sum(N(1:i))) = linspace(0,0,N(i));
    else
        for j=sum(N(1:i-1))+1:sum(N(1:i))
            kx(j) = kx(j-1) + (Time(j)-Time(j-1))*gamma/(2*pi)*Gx(i);
            ky(j) = ky(j-1) + (Time(j)-Time(j-1))*gamma/(2*pi)*Gy(i);
        end
    end
end




end

