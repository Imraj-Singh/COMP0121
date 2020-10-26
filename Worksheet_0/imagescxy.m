function h = imagescxy(C)
% function h = imagescxy(C)
%
% This function is a variation from imagesc(...) to make it behave
% differently in a number of ways:
%
% 1. always create a new figure
%
% 2. use the Cartesian axis mode, instead of the "matrix" axis mode as in
%    imagesc
%
% 3. use the aspect ratio that respects dimension
%
% 4. use the grayscale color map
%
% 4. include the colorbar
%
% author: Gary Hui Zhang (gary.zhang@ucl.ac.uk)
%

% create a new figure
h = figure;

% use imagesc to render the input image
imagesc(C);

% use the Cartesian axis mode; the system default is "axis ij"
axis xy;

% use the aspect ratio that respects dimension; this is equivalent to "axis
% equal" but additionally making the image use less space
axis image;

% use the grayscalar color map
colormap('gray');

% turn on the colorbar
colorbar;

end
