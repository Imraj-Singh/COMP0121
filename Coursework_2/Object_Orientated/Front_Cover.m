clc
clear
close all
z=5.5
x = linspace(0,z);

x1 = linspace(z,18);
y = sin(pi*x)./(pi*x);

y1 = sin(pi*(x1-2*z))./(pi*(x1-2*z));
plot([x x1],[y y1],'linewidth',2,'LineStyle','-','Color','k')
ylim([-1 1.5])
set(gca, 'Visible', 'off')


cmap = colormap(hsv(360));
alpha = atan2(W(j), U(j));
