% This script plots the radio stripes scenario. The APs are equally spaced
% in a circle, while the the UEs are contained inside an smaller one.
%
% This is version 1.0 (Last edited: 2025-04-29)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% monograph as described above.

%% SET PARAMETERS
rAPs = 60;      % Radius for APs
rUEs_max = 50;  % Limit radious for UEs
L = 48;         % Number of APs
K = 12;          % Number of UEs

%% GENERATE POSITIONS

%Positions of APs in a circle of radio r_APs equally spaced
theta_APs = linspace(0, 2*pi, L+1)';
cart_APs = rAPs .* (cos(theta_APs) + sin(theta_APs) * 1i);

%Random UEs locations with uniform distribution inside radio r_UEs
rUEs = rand(1, K) * rUEs_max;
theta_UEs =  rand(1, K) * 2 * pi;
cart_UEs = rUEs .* (cos(theta_UEs) + sin(theta_UEs) * 1i);

%% PLOT

% Extract real and imaginary parts
x_APs = real(cart_APs);
y_APs = imag(cart_APs);
x_UEs = real(cart_UEs);
y_UEs = imag(cart_UEs);

% Create the plot
figure;
hold on;
plot(x_APs, y_APs, 'ro', 'MarkerFaceColor', 'r'); % Red filled circles for APs
plot(x_UEs, y_UEs, 'ko', 'MarkerFaceColor', 'k'); % Black filled circles for UEs

% Draw circles
theta = linspace(0, 2*pi, 100);
x_circle_APs = rAPs * cos(theta);
y_circle_APs = rAPs * sin(theta);
x_circle_UEs = rUEs_max * cos(theta);
y_circle_UEs = rUEs_max * sin(theta);

plot(x_circle_APs, y_circle_APs, 'r--'); % Red dotted circle for APs
plot(x_circle_UEs, y_circle_UEs, 'k--'); % Black dotted circle for UEs

% Set axis equal for proper scaling
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
xlim([-rAPs-10, rAPs+10])
ylim([-rAPs-10, rAPs+10])

legend('APs', 'UEs', 'Location', 'northwest');
hold off;