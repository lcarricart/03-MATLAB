Ib = [1.04, 1.46, 2.55, 3.13, 3.99, 9.04, 12.82, 18.84];    % Base current I_B (µA)
Ic = [45,   59,   89,   103,   122,  209,   365,   446];   % Collector current I_C (µA)

figure
plot(Ib, Ic, '.-', 'LineWidth', 1, 'MarkerSize', 8)
xlabel('Base Current I_B (µA)')
ylabel('Collector Current I_C (µA)')
grid on
