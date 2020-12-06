% Small script to create front cover of coursework 2 - Comp MRI
% Imraj Singh 06/12/2020

clc
clear
close all
z=10.5;
x = linspace(0,z);

x1 = linspace(z,3*z);
y = sin(pi*x)./(pi*x);

y1 = sin(pi*(x1-2*z))./(pi*(x1-2*z));
plot([x x1],[y y1],'linewidth',2,'LineStyle','-','Color','k')
ylim([-1 1.5])
set(gca, 'Visible', 'off')
