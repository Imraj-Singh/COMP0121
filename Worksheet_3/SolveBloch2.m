function SolveBloch2(M0,omega0,T2)

% Define time span
tspan = [0 1000];
T1 = T2*10;
M01 = 1;
% Call ode45
[t,M] = ode45(@(t,M) [(M01-M(1))/T1; omega0 * M(3) - M(2)/T2; - omega0 * M(2) - M(3)/T2], tspan, M0);

% Plot
plot3(M(:,1),M(:,2),M(:,3),'-')
set(gca,'FontSize',30)
grid on
grid minor


end

