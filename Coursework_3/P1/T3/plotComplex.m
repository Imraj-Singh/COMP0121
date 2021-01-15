function [ h ] = plotComplex(xValues, yValues, subplotHandle)
% function [ h ] = plotComplex(xValues, yValues, subplotHandle)
%
% Plot complex-valued functions
%
% Input:
%
% xValues: the x coordinates of the points to plot
%
% yValues: the corresponding values of the function
%
% subplotHandle: an optional subplot handle
%
%
% OUTPUT:
%
% h - a subplot handle
%
%
% author: Dr Gary Zhang (gary.zhang@ucl.ac.uk)
%

% if a subplot handle is provided, using it; otherwise, create one
if exist('subplotHandle', 'var') ~= 0
    h = subplotHandle;
else
    h = subplot(1,1,1);
end

% plot the real values
plot(h, xValues, real(yValues), 'b.-');
% overlay the imagingary values
hold(h, 'on');
plot(h, xValues, imag(yValues), 'r.-');
% include the legend
legend(h, 'real', 'imag');

end

