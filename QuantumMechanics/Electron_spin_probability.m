% Very simple simulation which shows the probability of measuring the
% spin of an electron to be up realtive to a measurement device. The device is rotated
% over 360 degrees (from the inital starting position).
clear; clc;

% Unit vector for the orientation of electron spin (3D space)
n = [0.7071; 0.0; 0.7071];

% Unit vector for the orientation of the measurement device (3D space)
m = [0.7071; 0.0; 0.7071];


% Array of angles from 0 to 2pi
angles = linspace(0, 2*pi, 1000);

% Calculate the probability for each angle
average = dot(n, m) * cos(angles);
probabilities = (1 + average) / 2;

% Plotting
figure;
subplot(2,1,1);
plot(angles, probabilities);
title('Probability of Measuring Spin Up vs. Angle');
xlabel('Angle (radians)');
ylabel('Probability');
ylim([0, 1]);

subplot(2,1,2);
plot(angles,average);
title('Average value (for repeated experiments) vs. angle');
xlabel('Angle (radians)');
ylabel('Average');

% Display the maximum probability and corresponding angle
[max_prob, max_idx] = max(probabilities);
disp(['Maximum Probability: ' num2str(max_prob)]);
disp(['Corresponding Angle: ' num2str(angles(max_idx))]);