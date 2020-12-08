% Small script to create colour map for coursework 2 - Comp MRI
% Imraj Singh 08/12/2020
clc
clear
close all

% Set range of angles of spins
alpha = linspace(0,pi*2,361);

% Create circular colour map
cmap = colormap(hsv(360));
for i=1:length(alpha)
    % Plot a quiver in each of the angles
    h = quiver(0,0,cos(alpha(i)),sin(alpha(i)));
    
    % Convert the angle to degrees and round
    idx = ceil(rad2deg(alpha(i)));
    
    % If the degree are negative/zero then invalid and add 360 degrees
    if idx < 1
        idx = idx + 360;
    end
    
    % Set the colour of the quiver according to the colour map
    set(h, 'Color', cmap(idx,:))
    hold on
    % Square aspect ratio
    pbaspect([1 1 1])
    % Axes limits
    xlim([-1 1]);
    ylim([-1 1]);
    % Axes labels
    xlabel("$x'$", "interpreter", "latex", "fontsize", 15)
    ylabel("$y'$", "interpreter", "latex", "fontsize", 15)
    % Axes ticks
    xticks([-1 -0.5 0 0.5 1])
    yticks([-1 -0.5 0 0.5 1])
end

