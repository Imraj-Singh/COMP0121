clc
clear

%% Free procession with relaxation
% Author: Imraj Singh 03/11/2020

% Given parameters

% Prescribe T1
T2 = 1;
T1 = T2 * 2;

% Precession frequency
Onum = 100;
tsamp = 1001;
omega0 = 50 * pi;
deltaomega = 0;%omega0/100;
delta = 10;
omegaovector = omega0+deltaomega+delta*tan(pi*(rand(Onum,1)-1/2));

time = linspace(0, .5, tsamp);
signal = zeros(tsamp,Onum);
signal2 = zeros(tsamp,Onum);
signalres = zeros(tsamp,1);
signalres2 = zeros(tsamp,1);

for i=1:tsamp
    signal(i,:) = exp(-time(i)/T2)*sin(omegaovector*time(i));
    signalres(i) = sum(signal(i,:));
end

%plot(time,signal/Onum)
time3 = linspace(.5, 3, (tsamp*3-2));
time2 = linspace(.5, -2, (tsamp*3-2));
for i=1:(tsamp*3-2)
    signal2(i,:) = exp(-time3(i)/T2)*sin(omegaovector*time2(i));
    signalres2(i) = sum(signal2(i,:));
end

FinalSignal = [signalres; signalres2];
FinalTime = [time, time3];

figure
plot(FinalTime,FinalSignal/Onum)
hold on
plot(FinalTime,exp(-FinalTime/T2))



