% Small script to create front cover of coursework 2 - Comp MRI
% Imraj Singh 06/12/2020
clc
clear
close all

% This in an illustrative example of gradient echo the actual values mean
% nothing

% Variable to change to make the graph look nice
z=10.5;

% Set the x domain of initial sinc and calculate
x = linspace(0,z);
y = sin(pi*x)./(pi*x);

% Set x domain of second sinc and calculate
x1 = linspace(z,3*z);
y1 = sin(pi*(x1-2*z))./(pi*(x1-2*z));

% Plot the values
plot([x x1],[y y1],'linewidth',2,'LineStyle','-','Color','k')
% Set the y limits
ylim([-1 1.5])
% Remove the axes so that it is purely illustrative
set(gca, 'Visible', 'off')
