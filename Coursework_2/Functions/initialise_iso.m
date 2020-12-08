function iso = initialise_iso(x,y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
iso(1:length(x)) = struct();
for i=1:length(x)
        iso(i).x = x(i);
        iso(i).y = y(i);
        iso(i).Mx = [];
        iso(i).My = [];
        iso(i).Mz = [];
end
end

