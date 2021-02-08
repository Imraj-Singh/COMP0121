function iso = initialise_iso(x,y)
%INITALISE_ISO by Imraj Singh
% Inputs
%   x,y
% Outputs
%   iso

% Initialise the isochromat structure and set the x,y values 
iso(1:length(x)) = struct();
for i=1:length(x)
        iso(i).x = x(i);
        iso(i).y = y(i);
        iso(i).Mx = [];
        iso(i).My = [];
        iso(i).Mz = [];
        iso(i).Mxrec = [];
        iso(i).Myrec = [];
        iso(i).Mzrec = [];
end
end

