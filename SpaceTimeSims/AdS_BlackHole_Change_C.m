% Simulation itterates over values of C, creating a video of the different
% geometries

clear; clc;

%----------------------------------------
% PARAMETERS
%----------------------------------------
k = 1;  % curvature of AdS spacetime (positive value)
d = 3;  % Dimensionality of spacetime
C_values = 1:0.05:10;   % Define the range of C values

%----------------------------------------
% RECORDING
%----------------------------------------
% Set up the video writer
outputVideo = VideoWriter('AdS_Geometry_Visualization.mp4');
open(outputVideo);

for C = C_values
    r = (0:0.004:1)*2; % Radial coordinate
    t = (0:0.004:1)*2; % Adjusted omega range
    o = (0:0.004:1)*2*pi;

    % differentials
    dt = t(2) - t(1); % Time step
    dr = r(2) - r(1); % Radial step
    do = o(2) - o(1); % Angular step

    % Calculating Spacetime interval
    ds_squared = -(k^2 .* r.^2 + 1 - (C ./(r.^(d-2)))) .* dt^2 ...
        + 1 ./ (k^2 * r.^2 + 1 - (C ./ r.^(d-2))) .* dr^2 ...
        + r.^2 .* do^2;
    
    % Converting radial coordinates into Cartesian
    [R, Th] = meshgrid(r,o);
    x = R.*cos(Th);
    y = R.*sin(Th);
    [ds_squared] = meshgrid(ds_squared);
    
    % Create the 3D plot
    s = surf(x, y, ds_squared,'FaceAlpha',0.5);
    
    % Viewing position
    view(20, 25); % Camera tilt
    axis([-2,2, -2,2, -5e-3,20e-3]);
    
    title('Visualisation of AdS Black Hole');
    xlabel("x (spatial dimension)");
    ylabel("y (spatial dimension)");
    zlabel("Spacetime Interval (ds^2)");
    
    % Capture the frame and write to video
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
    
    % Clear the current figure for the next iteration
    clf;
end

% Close the video writer
close(outputVideo);
