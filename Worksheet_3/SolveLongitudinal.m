function SolveLongitudinal(Mz0,M0,T1)

% Define time span
tspan = [0 5];

% Call ode45
[t,Mz] = ode45(@(t,Mz) (M0-Mz)/T1, tspan, Mz0);

% Plot
plot(t,Mz,'-o')
xlabel("Time", "interpreter", "latex", "fontsize", 30)
ylabel("Magnetisation", "interpreter", "latex", "fontsize", 30)
set(gca,'FontSize',30)
grid on
grid minor

end

