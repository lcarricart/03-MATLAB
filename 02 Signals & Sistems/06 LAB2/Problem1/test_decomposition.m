% Script to test signal decomposition into even and odd parts
% Tests with the signal from Problem 1(a)

clear all;
close all;
clc;

% Define the signal from problem 1(a):
% x(t) = 0 for t < 0, x(t) = 0.5 for t = 0, x(t) = 1 for t > 0

% Create time vector (symmetric around 0)
N = 11; % Total number of points (must be odd)
n = -(N-1)/2:(N-1)/2; % Time indices from -(N-1)/2 to (N-1)/2

% Create signal vector
x = zeros(1, N);
center_idx = (N+1)/2; % Index of t=0

% Assign values according to problem definition
x(1:center_idx-1) = 0;     % t < 0
x(center_idx) = 0.5;        % t = 0
x(center_idx+1:end) = 1;    % t > 0

% Decompose the signal
[x_e, x_o] = decompose_signal(x);

% Verify: x = x_e + x_o
x_sum = x_e + x_o;

% Create plots
figure('Name', 'Signal Decomposition', 'Position', [100, 100, 1000, 800]);

% Plot 1: Original signal x(t)
subplot(2, 2, 1);
% Plot step function: 0 for t<0, 1 for t>0
hold on;
plot([min(n), 0], [0, 0], 'b', 'LineWidth', 2); % left part (t<0)
plot([0, max(n)], [1, 1], 'b', 'LineWidth', 2); % right part (t>0)
plot(0, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 1.5); % open circle at (0,0)
plot(0, 1, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 1.5); % open circle at (0,1)
plot(0, 0.5, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % filled circle at (0,0.5)
yline(0, 'k-', 'LineWidth', 0.5); % x-axis
xline(0, 'k-', 'LineWidth', 0.5); % y-axis
hold off;
grid on;
xlabel('Time (t)', 'FontSize', 11);
ylabel('Amplitude x(t)', 'FontSize', 11);
title('Original Signal x(t)');
axis([min(n)-0.5, max(n)+0.5, -0.2, 1.2]);
set(gca, 'XTick', n);

% Plot 2: Even part x_e(t)
subplot(2, 2, 2);
% Even part is constant at 0.5
plot([min(n), max(n)], [0.5, 0.5], 'r', 'LineWidth', 2);
hold on;
yline(0, 'k-', 'LineWidth', 0.5); % x-axis
xline(0, 'k-', 'LineWidth', 0.5); % y-axis
hold off;
grid on;
xlabel('Time (t)', 'FontSize', 11);
ylabel('Amplitude x_e(t)', 'FontSize', 11);
title('Even Part x_e(t) = 0.5[x(t) + x(-t)]');
axis([min(n)-0.5, max(n)+0.5, -0.2, 1.2]);
set(gca, 'XTick', n);

% Plot 3: Odd part x_o(t)
subplot(2, 2, 3);
% Odd part: -0.5 for t<0, 0 at t=0, 0.5 for t>0
hold on;
plot([min(n), 0], [-0.5, -0.5], 'g', 'LineWidth', 2); % left part (t<0)
plot([0, max(n)], [0.5, 0.5], 'g', 'LineWidth', 2); % right part (t>0)
plot(0, -0.5, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 1.5); % open circle at (0,-0.5)
plot(0, 0.5, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 1.5); % open circle at (0,0.5)
plot(0, 0, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % filled circle at (0,0)
yline(0, 'k-', 'LineWidth', 0.5); % x-axis
xline(0, 'k-', 'LineWidth', 0.5); % y-axis
hold off;
grid on;
xlabel('Time (t)', 'FontSize', 11);
ylabel('Amplitude x_o(t)', 'FontSize', 11);
title('Odd Part x_o(t) = 0.5[x(t) - x(-t)]');
axis([min(n)-0.5, max(n)+0.5, -0.6, 0.6]);
set(gca, 'XTick', n);

% Plot 4: Sum x_e(t) + x_o(t)
subplot(2, 2, 4);
% Sum should equal original: 0 for t<0, 1 for t>0, 0.5 at t=0
hold on;
plot([min(n), 0], [0, 0], 'm', 'LineWidth', 2); % left part (t<0)
plot([0, max(n)], [1, 1], 'm', 'LineWidth', 2); % right part (t>0)
plot(0, 0, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 1.5); % open circle at (0,0)
plot(0, 1, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 1.5); % open circle at (0,1)
plot(0, 0.5, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm'); % filled circle at (0,0.5)
yline(0, 'k-', 'LineWidth', 0.5); % x-axis
xline(0, 'k-', 'LineWidth', 0.5); % y-axis
hold off;
grid on;
xlabel('Time (t)', 'FontSize', 11);
ylabel('Amplitude x_e(t) + x_o(t)', 'FontSize', 11);
title('Sum: x_e(t) + x_o(t)');
axis([min(n)-0.5, max(n)+0.5, -0.2, 1.2]);
set(gca, 'XTick', n);

% Display results in command window
fprintf('\n=== Signal Decomposition Results ===\n\n');
fprintf('Time indices (n): ');
fprintf('%d ', n);
fprintf('\n\n');

fprintf('Original signal x(t): ');
fprintf('%.2f ', x);
fprintf('\n');

fprintf('Even part x_e(t):     ');
fprintf('%.2f ', x_e);
fprintf('\n');

fprintf('Odd part x_o(t):      ');
fprintf('%.2f ', x_o);
fprintf('\n');

fprintf('Sum x_e + x_o:        ');
fprintf('%.2f ', x_sum);
fprintf('\n\n');

% Verify decomposition
error = max(abs(x - x_sum));
fprintf('Maximum error |x - (x_e + x_o)|: %.10f\n', error);

if error < 1e-10
    fprintf('✓ Decomposition verified successfully!\n\n');
else
    fprintf('✗ Warning: Decomposition error is too large!\n\n');
end
