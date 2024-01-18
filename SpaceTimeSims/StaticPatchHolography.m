% Simulation of a black hole within de Sitter Space. Based on the space
% time metric used in static patch holography.
clear; clc;

%----------------------------------------
% PARAMETERS
%----------------------------------------
simRes = 0.005; % Defines the minimum distance between coordinates
M_sol = 100; % Mass of Black Hole (Solar Masses)

%----------------------------------------
% CALCULATIONS
%----------------------------------------
[rs, M] = calcBHRadius(M_sol);  % Finding Schwarzschild radius
R = 10*rs;  % Radius of the de Sitter space relative to black hole

% Generating radial coordiantes
[r, t, theta, phi, dr, dt, dtheta, dphi] = generateCoordinates(simRes, R); 

% Calculating spacetime interval
ds = spaceTimeInterval (r,R,M,theta,dt,dr,dtheta,dphi);

%----------------------------------------
% PLOT SPACETIME
%----------------------------------------
% Converting radial coordinates into Cartesian
[x, y] = sphToCart (r, theta);
[ds_matrix] = meshgrid(ds);

figure;
Spacetime = surf(x, y, ds_matrix, 'FaceAlpha', 0.9);

% Add title and labels
title('Visualization of Schwarzschild-de Sitter Black Hole');
xlabel('x (spatial dimension)');
ylabel('y (spatial dimension)');
zlabel('Spacetime Interval (ds)');

% Set viewing position
view(20, 25); % Camera tilt

%----------------------------------------
% PLOT HORIZONS (NOT FUNCTIONAL)
%----------------------------------------
% % Solve for roots of the horizon equation
% % horizon_equ = @(r) 1 - r.^3 / R^2 - 2*M*G;
% horizon_coefficients = [(-1/R^2), -2*M*G, 1];
% horizon_r = roots(horizon_coefficients);
% 
% % Filter real and positive roots
% horizon_r = horizon_r(imag(horizon_r) == 0 & real(horizon_r) > 0);
% 
% % Convert horizons to Cartesian coordinates
% x_horizon = horizon_r * sin(Theta_mesh);
% y_horizon = horizon_r * cos(Theta_mesh);
% z_horizon = horizon_r;
% 
% % Add horizons to the existing surf plot
% hold on;
% horizon = plot3(x_horizon, y_horizon, z_horizon);
% hold off;
% % Add legend for the horizon
% legend(horizon, 'Black Hole Horizon');

%----------------------------------------
% FUNCTIONS
%----------------------------------------
% Determines the radius of the event horizon according to the Schwarzschild Radius
function [rs, M] = calcBHRadius(M_sol)
    c = 2.99792458E8; % Speed of light
    G = 6.67408E-11; % Gravitational Constant
    M_sun = 1.98847e+30;  % mass of sun (kg)

    M = M_sol * M_sun; % Mass of Black Hole (Kg)
    rs = 2*G*M/c^2;  % Schwarzschild radius
end

% Generates spacetime coordinates and their intergal
function [r, t, theta, phi, dr, dt, dtheta, dphi] = generateCoordinates(simRes, R)

    % Create radial and angular coordinates
    r = (0:simRes:1) * R;
    t = (0:simRes:1) * R;
    theta = (0:simRes:1) * pi;
    phi = (0:simRes:1) * 2 * pi;

    % Calculate differential steps
    dr = diff(r);
    dt = diff(t);
    dtheta = diff(theta);
    dphi = diff(phi);

    % Duplicate the last element of each array to avoid incompatible sizes
    dr(end+1) = dr(end);
    dt(end+1) = dt(end);
    dtheta(end+1) = dtheta(end);
    dphi(end+1) = dphi(end);
end

% determines spacetime interval. (This is the degree of seperation between
% events within space time)
function ds = spaceTimeInterval (r,R,M,theta,dt,dr,dtheta,dphi)
    G = 6.67408E-11; % Gravitational Constant

    f_r = (1 - (r.^2 / R^2) - (2*M*G./r));
    ds_sqr = -f_r .* dt.^2 + (dr.^2 ./ f_r) .* r.^2 .* (dtheta.^2 + (sin(theta)).^2 .* dphi.^2);
    ds = sqrt(ds_sqr);
end

% Converts spherical coordinates to cartesian coordinates
function [x, y] = sphToCart (r, theta)
    [R_mesh, Theta_mesh] = meshgrid(r, theta);
    x = R_mesh .* sin(Theta_mesh);
    y = R_mesh .* cos(Theta_mesh);
end