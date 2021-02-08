function iso = initialise_B0(iso,B0)
%INITALISE_B0 by Imraj Singh
% Inputs
%   iso,B0
% Outputs
%   iso

% Set the B0 values
for i=1:length(iso(:,1))
    for j=1:length(iso(1,:))
        iso(i,j).B0 = B0;
    end
end
end

