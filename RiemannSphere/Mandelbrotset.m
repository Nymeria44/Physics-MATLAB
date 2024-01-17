clear; clc;

% Image size (pixels)
WIDTH = 600;
HEIGHT = 400;

% Plot window
RE_START = -2;
RE_END = 1;
IM_START = -1;
IM_END = 1;

% Create a grid of complex numbers corresponding to the Mandelbrot set
[X, Y] = meshgrid(linspace(RE_START, RE_END, WIDTH), linspace(IM_START, IM_END, HEIGHT));
Z = X + 1i*Y;

% Initialize arrays for Riemann Sphere and Cartesian coordinates
S_x = zeros(size(Z));
S_y = zeros(size(Z));
S_z = zeros(size(Z));

C_x = zeros(size(Z));
C_y = zeros(size(Z));

% Calculate the Mandelbrot set and map it to the Riemann Sphere
for i = 1:numel(Z)
    c = Z(i);
    z = 0;
    norm_z = 0;
    max_iter = 1000; % Maximum number of iterations (adjust as needed)

    % Iterate to determine if c is in the Mandelbrot set
    for iter = 1:max_iter
        z = z^2 + c;
        norm_z = abs(z);
        if norm_z > 2
            break; % Point is not in the Mandelbrot set
        end
    end
    
    % Check if the point is in the Mandelbrot set
    if norm_z <= 2
        % Calculate Riemann sphere coordinates
        S_x(i) = 2 * real(c) / (1 + abs(c)^2);
        S_y(i) = 2 * imag(c) / (1 + abs(c)^2);
        S_z(i) = (-1 + abs(c)^2) / (1 + abs(c)^2);

        % Calculate Cartesian coordinates
        C_x(i) = real(c);
        C_y(i) = imag(c);
    end
end

% Create a linear color map based on the value of c
c_values = linspace(1, 10, numel(Z));

% Create a 3D scatter plot for the Riemann Sphere coordinates
figure;

% Subplot for Riemann Sphere coordinates
subplot(1, 2, 1);
sphere % Generating reference sphere
hold on
% Plotting coordinates on the Riemann Sphere with colors
scatter3(S_x(:), S_y(:), S_z(:), [], c_values, 'filled');
% Set axis labels
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Riemann Sphere');
axis equal;
colorbar; % Add a colorbar to show the mapping

% Subplot for Cartesian coordinates
subplot(1, 2, 2);
hold on
% Plotting coordinates on the Cartesian grid
scatter(C_x(:), C_y(:), [], c_values, 'filled');
% Set axis labels
xlabel('Real');
ylabel('Imaginary');
title('Cartesian Grid');
axis equal;
colorbar; % Add a colorbar to show the mapping
