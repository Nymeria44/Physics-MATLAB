clear; clc;

% Define the parameters
c = 299792458; % Speed of light (m/s)
r_s = 0.2; % Schwarzschild radius
simRes = 0.005;

r = (0:simRes:1)*2; % Radial space
t = (0:simRes:1)*2; % time (reference frame of infinitely away from blackhole)
theta = (0:simRes:1)*pi; % Adjusted omega range
phi = (0:simRes:1)*2*pi;

% Differentals
dt = t(2) - t(1); % Time step
dr = r(2) - r(1); % Radial step
d_th = theta(2) - theta(1); % Angular step
d_phi = phi(2) - phi(1);

% Calculating Spacetime interval
g_o = d_th.^2 + (sin(theta)).^2 .* d_phi;
ds_sqr = - (1 - r_s./r) * c^2 .*dt.^2 + (1 - r_s./r).^(-1) .* dr.^2 + r.^2 .* g_o;


x = r*sin(theta)*cos(phi);
y = r.*sin(theta).*sin(phi);
z = r.*cos(theta);

% [x,y,z] = meshgrid(r.*sin(theta).*cos(phi), r.*sin(theta).*sin(phi), r.*cos(theta));
[ds_sqr] = meshgrid(ds_sqr);

% Create the 3D plot
Spacetime = surf(x,y,ds_sqr);

% Viewing position
% view(20, 25); % Camera tilt
% axis([-2,2, -2,2, -5e-3,20e-3]);

title('Visualisation of AdS Black Hole');
xlabel("x (spatial dimension)");
ylabel("y (spatial dimension)");
zlabel("Spacetime Interval (ds^2)");