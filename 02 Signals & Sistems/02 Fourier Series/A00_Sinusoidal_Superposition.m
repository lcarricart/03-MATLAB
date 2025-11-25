%% Purpose of the program: to perform wavefor synthesis of non-sinusoidal signals through sinusoidal superposition.
% This Fourier approximation leads to wavy values in the areas of interest (lets say I have squared signal, then instead of being straight, it's wavy). Gibb's over/under shootings are also present.

n = 99;         % Affects the amount of superimposed sine functions [Change this from 3 - 999]
f0 = 100;       % Frequency of the signal you want to resemble
t = [0:0.001:20]/1000;

%% Fourier descriptions of the functions Sawtooth and Triangle
for k = 1:2:n
    square(k,:)  = -4/(pi*k) * sin(2*pi*k*f0*t);
end

for i = 1:1:n
    sawtooth(i,:) = ((-2*(-1)^n) / (n*pi)) * sin(2*pi*i*f0*t);
    triangle(i,:) = 2/(pi^(2)*i^2) * sin(2*pi*i*f0*t);
end

% x = -pi:0.01:pi;
x = linspace(-2*pi, 2*pi, 1000);
function2 = @(x) sin(x) + cos(x);

%% Plotting
figure(1);
plot(t*1000, sum(square), 'r-', t*1000, sum(sawtooth), 'b-', t*1000, sum(triangle), 'g-', 'LineWidth', 2);
ylim([-2 2]);
title('Waveform Synthesis of Non-Sinusoidal Signals');
xlabel('t / ms');
ylabel('u / V');
title('Waveform Synthesis: Superposition of sinusoidal signals for N =', n);
grid on;
legend('Signal', 'Location', 'southeast');

%% Plotting of sin(x) + cos(x)
figure(2);
plot(x, function2(x));
title('Simple function addition plot');
grid on;