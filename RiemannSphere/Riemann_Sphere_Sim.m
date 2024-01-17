clear; clc;

% Create a grid of complex numbers covering the extended complex plane
grid = linspace(0, 30, 100000);

% [X, Y] = meshgrid(grid, grid);
X = grid;
Y = grid;
% Initialize arrays
S_x = zeros(size(X));
S_y = zeros(size(Y));
S_z = zeros(size(X));

C_x = zeros(size(X));
C_y = zeros(size(Y));

% Calculate the x, y, and z components for each point in the grid
for p = 1:numel(X)
    zeta = X(p) + i* sin(Y(p));
    norm_zeta = abs(zeta);
    
    % Point at infinity
    if norm_zeta == 0
        S_x(p) = 0;
        S_y(p) = 0;
        S_z(p) = -1;
        
        C_x(p) = 0;
        C_y(p) = 0;

    else
        % Calculating Riemann sphere
        S_x(p) = 2 * real(zeta) / (1 + zeta* conj(zeta));
        S_y(p) = 2 * imag(zeta) / (1 + zeta* conj(zeta));
        S_z(p) = (-1 + zeta * conj(zeta)) / (1 + zeta * conj(zeta));

        % Calculating Cartesian coordinates
        C_x(p) = real(zeta);
        C_y(p) = imag(zeta);
    end
end

% Create a 3D scatter plot for the Riemann sphere coordinates
figure;

% Subplot for Riemann sphere coordinates
subplot(1, 2, 1);
sphere % Generating reference sphere
hold on
% Plotting coordinates on the Riemann sphere
scatter3(S_x(:), S_y(:), S_z(:), 'filled');
% Set axis labels
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Riemann Sphere');
axis equal;

% Subplot for Cartesian coordinates
subplot(1, 2, 2);
hold on
% Plotting coordinates on the Cartesian grid
scatter(C_x(:), C_y(:), 'filled');
% Set axis labels
xlabel('Real');
ylabel('Img');
title('Cartesian Grid');
axis equal;
