clc
clear
close all
% author: Imraj Singh 06/01/2021

% Define the Field of view for the simulation
FOV = 8/1000;

% Calculate the k-space spacing
dk = 1/FOV;

% Define the number of points according to 2^p rule needed for IFFT
p = 7;
N = 2^p;

% Define the frequencies
k = linspace(-dk*(N/2),dk*(N/2-1),N);

% Define the distances
x = linspace(-FOV/2,FOV/2 - FOV/N,N);

% Set the video writing name and location
video = VideoWriter(['Ani', '.mp4'], 'MPEG-4');

% Set the frame rate of video
frameRate = 5;
video.set('FrameRate', frameRate);

% Open the video
video.open();
g = figure;

for i=1:N
    
    % Impulse location
    val = k(i);
    
    % Make signal vector with unit impulse
    signal = double(k == val);
    
    % original
    h = subplot(3,1,1);
    plotComplex(k/1000,signal,h);
    grid on
    box on
    hold off
    concate = append('Signal$ = \delta(k - $',int2str(round(k(i)/dk)),'$\Delta k)$');
    title(concate, "interpreter", "latex", "fontsize", 10)
    xlabel("$k$ (1/mm)", "interpreter", "latex", "fontsize", 10)
    ylabel("Signal", "interpreter", "latex", "fontsize", 10)
    
    % brute force
    h = subplot(3,1,2);
    if val == 0
        plotComplex(x*1000,exp(2*pi*1i*val*x),h);
    else
        plotComplex(x*1000,exp(2*pi*1i*nonzeros(signal.*k)*x),h);
    end
    grid on
    box on
    hold off
    title('Shift theorem', "interpreter", "latex", "fontsize", 10)
    % Set the axes labels
    xlabel("$x$ (mm)", "interpreter", "latex", "fontsize", 10)
    ylabel("$\rho$", "interpreter", "latex", "fontsize", 10)
    
    h = subplot(3,1,3);
    plotComplex(x*1000,fftshift(ifft(ifftshift(signal))*N),h);
    grid on
    box on
    hold off
    title('Brute force - IFFT(Signal) with shift', "interpreter", "latex", "fontsize", 10)
    % Set the axes labels
    xlabel("$x$ (mm)", "interpreter", "latex", "fontsize", 10)
    ylabel("$\rho$", "interpreter", "latex", "fontsize", 10)
    
    % Set frame to add to video
    frame = getframe(g);
    video.writeVideo(frame);
end

% Close the video
video.close();

% %% Brute force
% 
% out = zeros(1,N);
% 
% for ii=1:N
%     sum = 0.0;
%     for jj=1:N
%         sum = sum + kval(jj) * exp(2*1i*pi*k(ii)/dk*k(jj)/dk/N);
%     end
%     out(ii) = sum;
% end
%  
%  
% % shift of fourier of shift
% h = subplot(5,2,9:10);
% plotComplex(x,out,h);

