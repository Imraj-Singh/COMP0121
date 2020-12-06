function iso = initalise_B0(iso,B0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(iso(:,1))
    for j=1:length(iso(1,:))
        iso(i,j).B0 = B0;
    end
end
end

