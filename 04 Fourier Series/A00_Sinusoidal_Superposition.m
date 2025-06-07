%% Purpose of the program: to perform wavefor synthesis of non-sinusoidal signals through sinusoidal superposition.

n = 33;          % Affects the amount of superimposed sine functions [Change this from 3 - 999]
f0 = 100;       % Frequency of the signal you want to resemble
t = [0:0.001:20]/1000;

for k = 1:2:n
    square(k,:)  = -4/(pi*k) * sin(2*pi*k*f0*t);
end

%??????
for i = 1:1:n
    sawtooth(i,:) = ((-2*(-1)^n) / (n*pi)) * sin(2*pi*i*f0*t);
    triangle(i,:) = 2/(pi^(2)*i^2) * sin(2*pi*i*f0*t);
end

figure(1);
plot(t*1000, sum(square), 'r-', t*1000, sum(sawtooth), 'b-', t*1000, sum(triangle), 'g-', 'LineWidth', 2);
ylim([-2 2]);
title('Waveform Synthesis of Non-Sinusoidal Signals');
xlabel('t / ms');
ylabel('u / V');
title('Waveform Synthesis: Superposition of sinusoidal signals for N =', n);
grid on;
legend('Signal', 'Location', 'southeast');