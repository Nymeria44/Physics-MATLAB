% Simulation of a black hole in Anti-de Sitter space.

clear; clc;

%----------------------------------------
% PARAMETERS
%----------------------------------------
k = 0.5;  % Curvature of AdS spacetime (positive value)
C = 1;  % Describes the curvature of the black hole?
d = 3;  % Dimensionality of spacetime

simRes = 0.005; % Defines the minimum distance between coordinates

%----------------------------------------
% COORDINATES
%----------------------------------------
r = (0:simRes:1)*2; % Radial coordinate
t = (0:simRes:1)*2; % time (reference frame of infinitely away from blackhole)
o = (0:simRes:1)*2*pi; % Adjusted omega range

% Differentals
dt = t(2) - t(1); % Time step
dr = r(2) - r(1); % Radial step
do = o(2) - o(1); % Angular step

%----------------------------------------
% SPACETIME CALCULATION
%----------------------------------------
% Calculating Spacetime interval
ds_sqr = -(k^2 .* r.^2 + 1 - (C ./(r.^(d-2)))) .* dt^2 ...
    + 1 ./ (k^2 * r.^2 + 1 - (C ./ r.^(d-2))) .* dr^2 ...
    + r.^2 .* do^2;

%----------------------------------------
% PLOT
%----------------------------------------
% Converting radial coordinates into Cartesian
[R, Th] = meshgrid(r,o);
x = R.*cos(Th);
y = R.*sin(Th);
[ds_sqr] = meshgrid(ds_sqr);

% Create the 3D plot
Spacetime = surf(x, y, ds_sqr,'FaceAlpha',0.5);


% Viewing position
view(20, 25); % Camera tilt
axis([-2,2, -2,2, -5e-3,20e-3]);

title('Visualisation of AdS Black Hole');
xlabel("x (spatial dimension)");
ylabel("y (spatial dimension)");
zlabel("Spacetime Interval (ds^2)");