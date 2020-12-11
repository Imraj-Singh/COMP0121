function iso = initial_dir(iso,Axis,N)
%Initialise the spin direction
%   User defined direction of the spin, can only be set to
%   along the positive axis directions
if (nargin ~= 3)
    error('The number of arguments is: %d there should be only one argument either x, y or z.', nargin)
end
if any((Axis~='x')&(Axis~='y')&(Axis~='z'))
    error('Needs be an axis of either "x","y" or "z".')
end
if (Axis=='x')
    iso.Mx = [iso.Mx; ones(N,1)];
    iso.My = [iso.My; zeros(N,1)];
    iso.Mz = [iso.Mz; zeros(N,1)];
elseif (Axis=='y')
    iso.Mx = [iso.Mx; zeros(N,1)];
    iso.My = [iso.My; ones(N,1)];
    iso.Mz = [iso.Mz; zeros(N,1)];
elseif (Axis=='z')
    iso.Mx = [iso.Mx; zeros(N,1)];
    iso.My = [iso.My; zeros(N,1)];
    iso.Mz = [iso.Mz; ones(N,1)];
end
end

