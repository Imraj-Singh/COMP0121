% housekeeping
clc
clear

%% Rotation around the z-axis - visualisation
% Author - Imraj 11/11/2020

% initialise the video
video = VideoWriter(['2_1', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = 100;
video.set('FrameRate', frameRate);

% open video
video.open();

% initialise the figure
h = figure;

% define unit vectors
X = [1,0,0]';
Y = [0,1,0]';
Z = [0,0,1]';

% define angle of rotation 1 degree
theta = 2*pi/360;

% specify animation captures each degree of rotation
for i=0:360
    quiver3(0,0,0,X(1),X(2),X(3),'linewidth',2,'LineStyle','-','Color','k')
    hold on
    quiver3(0,0,0,Y(1),Y(2),Y(3),'linewidth',2,'LineStyle','--','Color','b')
    quiver3(0,0,0,Z(1),Z(2),Z(3),'linewidth',2,'LineStyle','-.','Color','r')
    
    % format title, legend, labels, limits, grid and box
    title(['z-axis rotation of unit vectors, $\theta =$', num2str(theta*i*360/(2*pi)),'$^\circ$'], "interpreter", "latex", "fontsize", 15)
    legend('$\hat{x}$','$\hat{y}$','$\hat{z}$', "interpreter", "latex", "fontsize", 15)
    xlabel("$M_x$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_y$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_z$", "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([-1 1]);
    hold off
    
    % assign frame and write it to the video
    frame = getframe(h);
    video.writeVideo(frame);
    
    % rotate by theta
    X = rotateZ(X,theta);
    Y = rotateZ(Y,theta);
    Z = rotateZ(Z,theta);
end

video.close();

%% Rotation around the x-axis - visualisation
% Author - Imraj

% initialise the video
video = VideoWriter(['2_2', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = 100;
video.set('FrameRate', frameRate);

% open video
video.open();

% initialise the figure
h = figure;

% define unit vectors
X = [1,0,0]';
Y = [0,1,0]';
Z = [0,0,1]';

% define angle of rotation 1 degree
theta = 2*pi/360;

% specify animation captures each degree of rotation
for i=0:360
    quiver3(0,0,0,X(1),X(2),X(3),'linewidth',2,'LineStyle','-','Color','k')
    hold on
    quiver3(0,0,0,Y(1),Y(2),Y(3),'linewidth',2,'LineStyle','--','Color','b')
    quiver3(0,0,0,Z(1),Z(2),Z(3),'linewidth',2,'LineStyle','-.','Color','r')
    
    % format title, legend, labels, limits, grid and box
    title(['x-axis rotation of unit vectors, $\theta =$', num2str(theta*i*360/(2*pi)),'$^\circ$'], "interpreter", "latex", "fontsize", 15)
    legend('$\hat{x}$','$\hat{y}$','$\hat{z}$', "interpreter", "latex", "fontsize", 15)
    xlabel("$M_x$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_y$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_z$", "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([-1 1]);
    hold off
    
    % assign frame and write it to the video
    frame = getframe(h);
    video.writeVideo(frame);
    
    % rotate by theta
    X = rotateX(X,theta);
    Y = rotateX(Y,theta);
    Z = rotateX(Z,theta);
end

video.close();

%% Rotation around the y-axis - visualisation
% Author - Imraj

% initialise the video
video = VideoWriter(['2_3', '.mp4'], 'MPEG-4');

% set the frame rate
frameRate = 100;
video.set('FrameRate', frameRate);

% open video
video.open();

% initialise the figure
h = figure;

% define unit vectors
X = [1,0,0]';
Y = [0,1,0]';
Z = [0,0,1]';

% define angle of rotation 1 degree
theta = 2*pi/360;

% specify animation captures each degree of rotation
for i=0:360
    quiver3(0,0,0,X(1),X(2),X(3),'linewidth',2,'LineStyle','-','Color','k')
    hold on
    quiver3(0,0,0,Y(1),Y(2),Y(3),'linewidth',2,'LineStyle','--','Color','b')
    quiver3(0,0,0,Z(1),Z(2),Z(3),'linewidth',2,'LineStyle','-.','Color','r')
    
    % format title, legend, labels, limits, grid and box
    title(['y-axis rotation of unit vectors, $\theta =$', num2str(theta*i*360/(2*pi)),'$^\circ$'], "interpreter", "latex", "fontsize", 15)
    legend('$\hat{x}$','$\hat{y}$','$\hat{z}$', "interpreter", "latex", "fontsize", 15)
    xlabel("$M_x$", "interpreter", "latex", "fontsize", 15)
    ylabel("$M_y$", "interpreter", "latex", "fontsize", 15)
    zlabel("$M_z$", "interpreter", "latex", "fontsize", 15)
    grid on
    box on
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([-1 1]);
    hold off
    
    % assign frame and write it to the video
    frame = getframe(h);
    video.writeVideo(frame);
    
    % rotate by theta
    X = rotateY(X,theta);
    Y = rotateY(Y,theta);
    Z = rotateY(Z,theta);
end

video.close();
