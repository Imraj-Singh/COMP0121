clc
clear

%% Free procession with relaxation
% Author: Imraj Singh 03/11/2020

% Given parameters

% Prescribe T1
T2 = 1;
T1 = T2 * 2;
delta = 50 * pi/100*2;

Onum = 100;
tsamp = 1001;
omega0 = 50 * pi;
time = linspace(0, .5, tsamp);
signal = zeros(tsamp,Onum);
signal2 = zeros(tsamp,Onum);
signalres = zeros(tsamp,1);
signalres2 = zeros(tsamp,1);
time2 = linspace(.5, -2, (tsamp*3-1));
time3 = linspace(.5, 3, (tsamp*3-1));

    

video = VideoWriter(['5_2sigDo', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = 0.25;
video.set('FrameRate', frameRate,'Quality',100);

video.open();

h = figure;



for z=[-omega0*.99 -omega0*.75 -omega0*.5 -omega0*.25 0 omega0*.25 omega0*.5 omega0 omega0*2]
    deltaomega = z;
    omegaovector = omega0+deltaomega+delta*tan(pi*(rand(Onum,1)-1/2));
    
    for i=1:tsamp
        signal(i,:) = exp(-time(i)/T2)*sin(omegaovector*time(i));
        signalres(i) = sum(signal(i,:));
    end
    
    for i=1:(tsamp*3-1)
        signal2(i,:) = exp(-time3(i)/T2)*sin(omegaovector*time2(i));
        signalres2(i) = sum(signal2(i,:));
    end
    
    FinalSignal = [signalres; signalres2];
    FinalTime = [time, time3];
    
    plot(FinalTime,FinalSignal/Onum)
    hold on
    plot(FinalTime,exp(-FinalTime/T2))
    title(['$\omega_0 = $', num2str(deltaomega/omega0),'$\omega$'], "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([0 3]);
    ylim([-1 1]);
    xlabel("Time (s)", "interpreter", "latex", "fontsize", 15)
    ylabel("Normalised Signal", "interpreter", "latex", "fontsize", 15)
    hold off
    frame = getframe(h);
    video.writeVideo(frame);
    
end

video.close();

