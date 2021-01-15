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

% Define Gaussian
gaus = @(k,mu,sig)(1/(sig*sqrt(2)))*exp(-(((k-mu).^2)/(2*sig.^2)));

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
    signal = gaus(k,0,8*dk).*double(k == val);
    % original
    h = subplot(3,1,1);
    plotComplex(k/1000,signal,h);
    hold on
    plot(k/1000, gaus(k,0,8*dk), '--k')
    xlim([k(1)/1000 k(end)/1000])
    ylim([0 max(gaus(x,0,8*dk))])
    grid on
    box on
    hold off
    concate = append('Signal$ = G(k)*\delta(k - $',int2str(round(k(i)/dk)),'$\Delta k)$');
    title(concate, "interpreter", "latex", "fontsize", 10)
    legend(h, 'real', 'imag', 'Gaussian');
    xlabel("$k$ (1/mm)", "interpreter", "latex", "fontsize", 10)
    ylabel("Signal", "interpreter", "latex", "fontsize", 10)
    
    % brute force
    h = subplot(3,1,2);
    if val == 0
        plotComplex(x*1000,nonzeros(signal)*exp(2*pi*1i*x),h);
    else
        plotComplex(x*1000,nonzeros(signal)*exp(2*pi*1i*nonzeros(k.*double(k == val))*x),h);
    end
    grid on
    box on
    hold off
    xlim([x(1)*1000 x(end)*1000])
    lims_val = fftshift(ifft(ifftshift(gaus(k,0,8*dk).*[zeros(1,64),1,zeros(1,63)]))*N);
    ylim([-max(lims_val) max(lims_val)])
    title('Linearity', "interpreter", "latex", "fontsize", 10)
    % Set the axes labels
    xlabel("$x$ (mm)", "interpreter", "latex", "fontsize", 10)
    ylabel("$\rho$", "interpreter", "latex", "fontsize", 10)
    
    h = subplot(3,1,3);
    plotComplex(x*1000,fftshift(ifft(ifftshift(signal))*N),h);
    grid on
    box on
    hold off
    xlim([x(1)*1000 x(end)*1000])
    lims_val = fftshift(ifft(ifftshift(gaus(k,0,8*dk).*[zeros(1,64),1,zeros(1,63)]))*N);
    ylim([-max(lims_val) max(lims_val)])
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

