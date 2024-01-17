% Clear workspace and command window
clear; clc;

%----------------------------------------
% PARAMETERS
%----------------------------------------
% Constants
G = 6.67408E-11; % Gravitational Constant
c = 2.99792458E8; % Speed of light
M_sun = 1.98847e+30;  % mass of sun (kg)

% Defining parameters
simRes = 0.005; % Resolution of the coordinates

M_sol = 0.1; % Mass of Black Hole (Solar Masses)
M = M_sol * M_sun; % Mass of Black Hole (Kg)
rs = 2*G*M/c^2;  % Schwarzschild radius
R = 10*rs;  % Radius of the de Sitter space relative to black hole

%----------------------------------------
% COORDINATES
%----------------------------------------
% Creating radial and angular coordinates
r = (0:simRes:1)*R;
t = (0:simRes:1)*R;
theta = (0:simRes:1)*pi;
phi = (0:simRes:1)*2*pi;

% Calculate differential steps
dr = diff(r);
dt = diff(t);
dtheta = diff(theta);
dphi = diff(phi);

% Duplicating last element of array to avoid incompatible sizes
dr(end+1) = r(end) - r(end-1);
dt(end+1) = t(end) - t(end-1);
dtheta(end+1) = theta(end) - theta(end-1);
dphi(end+1) = phi(end) - phi(end-1);

%----------------------------------------
% SPACETIME CALC
%----------------------------------------
% Calculate spacetime interval
f_r = (1 - (r.^2 / R^2) - (2*M*G./r));
ds_sqr = -f_r .* dt.^2 + (dr.^2 ./ f_r) .* r.^2 .* (dtheta.^2 + (sin(theta)).^2 .* dphi.^2);
ds = sqrt(ds_sqr);

% Convert radial coordinates into Cartesian
[R_mesh, Theta_mesh] = meshgrid(r, theta);
x = R_mesh .* sin(Theta_mesh);
y = R_mesh .* cos(Theta_mesh);
[ds_matrix] = meshgrid(ds);

%----------------------------------------
% PLOT
%----------------------------------------
figure;
Spacetime = surf(x, y, ds_matrix, 'FaceAlpha', 0.5);

% Add title and labels
title('Visualization of Schwarzschild-de Sitter Black Hole');
xlabel('x (spatial dimension)');
ylabel('y (spatial dimension)');
zlabel('Spacetime Interval (ds^2)');

% Identify and plot horizons
% horizon_r = r(abs(f_r) == min(abs(f_r))); % Find the horizon location
% horizon = plot3(horizon_r * sin(pi), horizon_r * cos(pi), 0, 'ro', 'MarkerSize', 10);

% Set viewing position and axis limits
view(20, 25); % Camera tilt


% Add legend for the horizon
% legend(horizon, 'Black Hole Horizon');