% MATLAB Script for Poincaré Disk Model finding geodesic
clear; clc;

%----------------------------------------
% PARAMETERS
%----------------------------------------
boundaryRes = 100; % Number of points used to generate boundary of disk
A = [0.1, 0.1]; % Point within the unit circle
O = [0, 0]; % Origin

%----------------------------------------
% Calculations
%----------------------------------------
% Creating boundary of the disk
Gamma = createCircle(O, 1, boundaryRes);

% Creating symbolic variables


% Creating line passing through O A
OA_line = createLineSymbolic(O, A);

% Finding the inversion point B
B = findInversionPoint(A,O);

%----------------------------------------
% Plotting
%----------------------------------------
% Plotting the boundary circle
hold on;
plot(Gamma(:, 1), Gamma(:, 2), 'k-', 'DisplayName', 'Unit Circle');

% Plotting the line OA
fplot(OA_line, 'b--', 'DisplayName', 'Line OA');

% Plotting the point A
plot(A(1), A(2), 'ro', 'DisplayName', 'Point A');

% Plotting the inversion point B
plot(B(1), B(2), 'go', 'DisplayName', 'Point B');
grid on;
xlim([-5 5])
ylim([-5 5])

%----------------------------------------
% FUNCTIONS
%----------------------------------------
% Creates boundary of the disk (Unit circle)
function circle = createCircle(O, r, boundaryRes)
    theta = linspace(0, 2*pi, boundaryRes);
    circle = O + r * [cos(theta); sin(theta)]';
end

% Creates a symbolic line passing through points O and A
function line = createLineSymbolic(O, A)
    syms x;    
    % Check if the line is vertical
    if (A(1) - O(1)) == 0
        line = x - O(1);
    else
        % Calculate the symbolic equation for the line passing through O and A
        m = (A(2) - O(2)) / (A(1) - O(1)); % slope
        b = O(2) - m * O(1); % y-intercept
        line = m * x + b;
    end
end

% Finds the point of inversion
function B = findInversionPoint(A,O)
    % Calculate length of OB
    mag_OB = 1 / norm(A);

    % Determine the direction of OA
    direction_OA = (A - O) / norm(A - O);

    % Calculate the coordinates of B along the line OA
    B = mag_OB * direction_OA;
end
