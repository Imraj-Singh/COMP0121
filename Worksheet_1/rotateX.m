function [rotated] = rotateZ(vector,theta)
%Rotate around z-axis
%   vector - the vector you want rotated, must have three dimensions
%   theta - angle in radians to rotate the vector around the z-axis
%   clockwise
%   

rotated(1) = vector(1)*cos(theta) + vector(2)*sin(theta);
rotated(2) = vector(2)*cos(theta) - vector(1)*sin(theta);
rotated(3) = vector(3);
end

