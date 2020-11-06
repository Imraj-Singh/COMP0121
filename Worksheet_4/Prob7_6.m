clc
clear

domega = 5;

omegao = 100;

time = linspace(0,3,10001);

modulated = sin((omegao+domega).*time) + sin((omegao-domega).*time);

demodulated1 = modulated.*sin((omegao).*time);
demodulated2 = modulated.*sin((omegao+domega).*time);
% demodulated2 = sin(omegao-domega).*time.*sin((omegao-domega).*time);

hold on
plot(time,modulated)
%plot(time,demodulated1)
plot(time,demodulated2)




