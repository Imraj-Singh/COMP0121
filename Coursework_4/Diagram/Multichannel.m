clc
clear

I = imread('Peppers.png');

[r,g,b] = imsplit(I);

montage({r,g,b},'Size',[1 3])