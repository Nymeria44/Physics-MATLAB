% --------------------------------------------------------------------------------
% MATLAB Script for finding geodesic on a Poincaré Disk Model (Hyberbolic
% Geometry). 
% 
% This is done by taking a point, A, finding the inverse of the
% point, B, (located outside of the unit sphere). Then the circle passing
% through both A and B will produce an arc inside of the Poincaré Disk, tracing 
% a Geodesic (straight line) in the Hyperbolic Geometry.
%--------------------------------------------------------------------------------
clear; clc;

%--------------------------------------------------------------------------------
% PARAMETERS
%--------------------------------------------------------------------------------
simRes = 100; % Number of points used to generate lines
debug = false;
A = [0.5, 0]; % Point within the unit circle
O = [0, 0]; % Origin

%--------------------------------------------------------------------------------
% CALCULATIONS
%--------------------------------------------------------------------------------
% Creating boundary of the Poincaré Disk
Gamma = createCircle(O, 1, simRes);

% Creating line passing through O A
OA_line = createSymbolicLine(O, A);

% Finding the inversion point B
B = findInversionPoint(A,O);

% Finding midpoint between A and B
AB_mid = (A + B) / 2;

% Finding radius for circle passing through A and B
AB_r = sqrt((A(1) - B(1))^2 + (A(2) - B(2))^2)/2;

% Creating circle which traces the geodesic
AB_circ = createCircle(AB_mid, AB_r, simRes);

% Determing the points where Poincaré Disk and AB circle intersect
[I1, I2] = findCircleIntersections(O, 1, AB_mid, AB_r);

% Finding Geodesic of disk
Geodesic = createGodesicArc (I1, I2, AB_mid, AB_r, simRes);

%--------------------------------------------------------------------------------
% PLOTTING
%--------------------------------------------------------------------------------
% Plotting the boundary of the Poincaré Disk
hold on;
plot(Gamma(:, 1), Gamma(:, 2), 'k-', 'DisplayName', 'Unit Circle');

% Plotting AB circle for debugging
plot(AB_circ(:, 1), AB_circ(:, 2), 'k-', 'DisplayName', 'AB circle');

% Plotting the geodesic line
plot(Geodesic(:, 1), Geodesic(:, 2), 'k-', 'DisplayName', 'Geodesic Line')

% Plotting the line OA
fplot(OA_line, 'b--', 'DisplayName', 'Line OA');

% Plotting the point A
plot(A(1), A(2), 'ro', 'DisplayName', 'Point A');
text(A(1), A(2), '  A', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Plotting the inversion point B
plot(B(1), B(2), 'go', 'DisplayName', 'Point B');
text(B(1), B(2), '  B', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% Plotting midpoint of AB
plot(AB_mid(1), AB_mid(2), 'co', 'DisplayName', 'Midpoint of AB');
text(AB_mid(1), AB_mid(2), '  Midpoint of AB', 'HorizontalAlignment', ...
    'left', 'VerticalAlignment', 'bottom');

% Plotting points of intersection between circles
plot(I1(1), I1(2), 'mo', 'DisplayName', 'Intersection P1');
text(I1(1), I1(2), '  Intersection P1', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'bottom');

plot(I2(1), I2(2), 'mo', 'DisplayName', 'Intersection P2');
text(I2(1), I2(2), '  Intersection P2', 'HorizontalAlignment', ...
    'left', 'VerticalAlignment', 'bottom');

grid on;
axis equal;
legend('Location', 'best');

%--------------------------------------------------------------------------------
% FUNCTIONS
%--------------------------------------------------------------------------------
% Creates boundary of the disk (Unit circle)
function circle = createCircle(O, r, boundaryRes)
    theta = linspace(0, 2*pi, boundaryRes);
    circle = O + r * [cos(theta); sin(theta)].';
end

% Creates a symbolic line passing through points O and A
function line = createSymbolicLine(O, A)
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

% Finds intersection between circles
function [intersection1, intersection2] = findCircleIntersections(O1, r1, O2, r2)
    % Symbolic variables
    syms x y;

    % Equations for the circles
    circle1 = (x - O1(1))^2 + (y - O1(2))^2 - r1^2;
    circle2 = (x - O2(1))^2 + (y - O2(2))^2 - r2^2;

    % Solve the system of equations
    sol = solve([circle1, circle2], [x, y]);

    % Extract intersection points from the solution
    intersection1 = double([sol.x(1) + O1(1), sol.y(1) + O1(2)]);
    intersection2 = double([sol.x(2) + O1(1), sol.y(2) + O1(2)]);
end

function arc = createGodesicArc(P1, P2, O, r, boundaryRes, debug)
    % Validate input parameters
    if nargin < 6
        debug = false;
    end
    
    if ~isnumeric(boundaryRes) || ~isscalar(boundaryRes) || boundaryRes <= 0
        error('boundaryRes must be a positive scalar.');
    end

    % Calculate angles
    angle_P1 = atan2d(P1(2) - O(2), P1(1) - O(1));
    angle_P2 = atan2d(P2(2) - O(2), P2(1) - O(1));

    % Handle angle wrapping
    angle_diff = angle_P2 - angle_P1;
    if angle_diff < -180
        angle_diff = angle_diff + 360;
    elseif angle_diff > 180
        angle_diff = angle_diff - 360;
    end

    % Generate parameter values
    t = linspace(0, angle_diff, boundaryRes);

    % Parametric equations for the arc using circr function
    arc = O + r * [cosd(angle_P1 + t); sind(angle_P1 + t)].';

    % Filter points outside the unit circle
    arc = arc(vecnorm(arc') <= 1, :);

    % Optionally display arc dimensions (for debugging)
    if debug
        disp('Arc dimensions:')
        disp(['arc size:', num2str(size(arc))])
    end
end