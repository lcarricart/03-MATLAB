x = [669 700 740 790 821 852];                  % Drain Source Voltage (mV)
y = [0.080 0.164 0.382 0.880 1.170 1.960];      % Drain Current (mA)

figure
plot(x, y, '.-', 'LineWidth', 4, 'MarkerSize', 8)
xlabel('Drain Source Voltage (mV)')
ylabel('Drain Current (mA)')
grid on
