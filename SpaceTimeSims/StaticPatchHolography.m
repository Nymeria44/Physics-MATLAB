% Simulation of a black hole within de Sitter Space. Based on the space
% time metric used in static patch holography.
clear; clc;

%----------------------------------------
% PARAMETERS
%----------------------------------------
simRes = 0.005; % Defines the minimum interval between coordinates (from 0 to 1)

% NOTE: the spacetime metric is only valid for small black holes
M = 4000; % % Mass of Black Hole (KG)

%----------------------------------------
% CALCULATIONS
%----------------------------------------
rs = calcBHRadius(M);  % Finding Schwarzschild radius
R = rs * 10E10;  % Radius of the de Sitter space relative to black hole

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
Spacetime = surf(x, y, ds_matrix);
colormap('jet');

% Dotted lines on mesh for better visibility
set(Spacetime,'LineStyle',':') 

%----------------------------------------
% PLOT HORIZONS
%----------------------------------------
% finding horizon of de Sitter space
r_h = findHorizons(R, M);
ds_h = spaceTimeInterval (r_h,R,M,theta,dt,dr,dtheta,dphi);
[x_h, y_h] = sphToCart (r_h, theta);

hold on;
deSitterHorizon = plot3(x_h, y_h, ds_h, 'r', 'LineWidth', 4);
hold off;

% legend, title and labels
legend(deSitterHorizon, 'de Sitter Horizon');
title('Visualization of Schwarzschild-de Sitter Black Hole');
xlabel('x (spatial dimension)');
ylabel('y (spatial dimension)');
zlabel('Spacetime Interval (ds)');

% Set viewing position
view(20, 25); % Camera tilt

%----------------------------------------
% FUNCTIONS
%----------------------------------------
% Determines the radius of the event horizon according to the Schwarzschild Radius
function [rs] = calcBHRadius(M)
    c = 2.99792458E8; % Speed of light
    G = 6.67408E-11; % Gravitational Constant

    rs = 2*G*M/c^2;  % Schwarzschild radius
end

% Generates spacetime coordinates and their intergal
function [r, t, theta, phi, dr, dt, dtheta, dphi] = generateCoordinates(simRes, R)
    % Create radial and angular coordinates
    r = (0:simRes:1) * R;
    t = (0:simRes:1) * (R/2);
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

% finds value of r at horizons
function r_h = findHorizons(R, M)
    G = 6.67408E-11; % Gravitational Constant

    % Solving f_r equation for roots
    horizon_coefficients = [(-1/R^2), -2*M*G, 1];
    r_h = roots(horizon_coefficients);

    % Filtering for real and positive roots
    r_h = r_h(imag(r_h) == 0 & real(r_h) > 0);
end