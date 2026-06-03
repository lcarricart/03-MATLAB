% PLOT_FOURIER_SERIES
% Script to compute and plot Fourier series of a square wave
% Problem 2(b): Plots the Fourier series reconstruction

clear all;
close all;
clc;

%% Get parameters from user via input dialogs
prompt = {'Number of Fourier coefficients (N):', ...
          'Fundamental frequency f0 [Hz]:', ...
          'Start time [s]:', ...
          'End time [s]:', ...
          'Number of time points:'};
dlgtitle = 'Fourier Series Parameters';
dims = [1 60];
defaults = {'10', '1', '-2', '2', '1000'};
answer = inputdlg(prompt, dlgtitle, dims, defaults);

% Check if user cancelled
if isempty(answer)
    disp('Operation cancelled by user.');
    return;
end

% Parse user input
N = str2double(answer{1});
f0 = str2double(answer{2});
t_start = str2double(answer{3});
t_end = str2double(answer{4});
num_points = str2double(answer{5});

% Validate input
if isnan(N) || isnan(f0) || N < 1 || f0 <= 0
    error('Invalid input! N must be >= 1 and f0 must be > 0.');
end

if isnan(t_start) || isnan(t_end) || isnan(num_points) || t_end <= t_start || num_points < 10
    error('Invalid time vector! t_end must be > t_start and num_points must be >= 10.');
end

N = round(N);           % Ensure N is an integer
num_points = round(num_points);  % Ensure num_points is an integer
T0 = 1/f0;              % Period [seconds]

% Create time vector based on user input
t = linspace(t_start, t_end, num_points);

%% Compute Fourier series
[a_k, b_k, x_fourier] = fourier_square_wave(N, f0, t);

%% Generate ideal square wave for comparison
x_ideal = square(2*pi*f0*t);  % MATLAB's square wave (-1 to +1)

%% Create plots
figure('Name', 'Fourier Series of Square Wave', 'Position', [100, 100, 1400, 600]);

% Plot 1: Original (Ideal) Square Wave
subplot(1, 2, 1);
hold on;
plot(t, x_ideal, 'b-', 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 0.5);
xline(0, 'k-', 'LineWidth', 0.5);
hold off;
grid on;
xlabel('Time (t) [s]', 'FontSize', 12);
ylabel('Amplitude x(t)', 'FontSize', 12);
title('Original Square Wave', 'FontSize', 13);
axis([min(t), max(t), -1.5, 1.5]);

% Plot 2: Fourier Series Approximation
subplot(1, 2, 2);
hold on;
stem(t, x_fourier, 'r', 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 0.5);
xline(0, 'k-', 'LineWidth', 0.5);
hold off;
grid on;
xlabel('Time (t) [s]', 'FontSize', 12);
ylabel('Amplitude x(t)', 'FontSize', 12);
title(sprintf('Fourier Series Approximation (N=%d, f_0=%.1f Hz)', N, f0), 'FontSize', 13);
axis([min(t), max(t), -1.5, 1.5]);

%% Display coefficients in command window
fprintf('\n=== Fourier Series Coefficients ===\n');
fprintf('Fundamental frequency f0 = %.2f Hz\n', f0);
fprintf('Period T0 = %.4f s\n', T0);
fprintf('Number of coefficients N = %d\n', N);
fprintf('Time vector: t = [%.2f, %.2f] s with %d points\n\n', t_start, t_end, num_points);

fprintf('DC component a_0 = %.4f\n\n', a_k(1));

fprintf('k\ta_k\t\tb_k\t\t|c_k|\n');
fprintf('-------------------------------------------\n');
for k = 1:N
    mag = sqrt(a_k(k+1)^2 + b_k(k)^2);
    fprintf('%d\t%.6f\t%.6f\t%.6f\n', k, a_k(k+1), b_k(k), mag);
end
fprintf('\n');
