clear; clc;
% Parameters
r = 0.5;   % Radius of the black hole
k = 0;   % Constant curvature of spacetime
gridDimension = -6:0.05:6;

% Creating 2D Spacetime
[x, y] = meshgrid(gridDimension, gridDimension);

% Calculate distance function (using inverse of distance)
d = k - 1 ./ sqrt(x.^2 + y.^2 - r^2);

% Create a 3D meshgrid for embedding
[X, Y] = meshgrid(gridDimension, gridDimension);
Z = sqrt(X.^2 + Y.^2); % Embedding in 3D space

% % Plot the visualization
figure;
% surf(X, Y, Z);
% hold on;
scatter3(x(:), y(:), d(:), 'filled'); % Scatter plot of points

% % Plot the event horizon (a cylinder)
% cylinder_radius = r;  % Radius of the cylinder representing the event horizon
% height = 3;             % Height of the cylinder in the z-direction
% theta = linspace(0, 2*pi, 100);
% z_event_horizon = linspace(0, height, 100);
% [X_cylinder, Z_cylinder] = meshgrid(cylinder_radius * cos(theta), z_event_horizon);
% Y_cylinder = cylinder_radius * sin(theta);
% surf(X_cylinder, Y_cylinder, Z_cylinder, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'FaceColor', 'blue');

xlabel('x');
ylabel('y');
zlabel('z');
title('Embedding of 2D Space with Spherical Void and Event Horizon into 3D Space');
