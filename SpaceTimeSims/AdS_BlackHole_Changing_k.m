% Simulation itterates over values of k, creating a video of the different
% geometries

clear; clc;

%----------------------------------------
% PARAMETERS
%----------------------------------------
C = 1;  % Describes the curvature of the black hole?
d = 3;  % Dimensionality of spacetime
K_values = [0.001:0.01:1, 1:1:100]; % Define the range of k values

%----------------------------------------
% RECORDING
%----------------------------------------
% Set up the video writer
outputVideo = VideoWriter('AdS_Geometry_Visualization.mp4');
open(outputVideo);

for k = K_values
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
    
    % Tilt the camera
    view(20, 20);
    
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
