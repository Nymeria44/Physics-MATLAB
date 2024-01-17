% Define the vectors n and m
n = [0.7071; 0.0; 0.7071]; % For example, n = [sqrt(2)/2; 0; sqrt(2)/2]
m = [0.5; 0.5; 0.7071];   % For example, m = [0.5; 0.5; sqrt(2)/2]

% Define an array of angles from 0 to pi
angles = linspace(0, pi, 1000);

% Initialise an array to store probabilities
probabilities = zeros(size(angles));

% Calculate the probability for each angle
for i = 1:length(angles)
    angle = angles(i);
    cos_theta_mn = dot(n, m) * cos(angle);
    probabilities(i) = (1 + cos_theta_mn) / 2;
end

% Create a plot
figure;
plot(angles, probabilities);
title('Probability of Measuring Spin Up vs. Angle');
xlabel('Angle (radians)');
ylabel('Probability');
ylim([0, 1]);

% Display the maximum probability and corresponding angle
[max_prob, max_idx] = max(probabilities);
disp(['Maximum Probability: ' num2str(max_prob)]);
disp(['Corresponding Angle: ' num2str(angles(max_idx))]);
