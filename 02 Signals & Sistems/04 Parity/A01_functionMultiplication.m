%% Purpose of the program: define even and odd functions to demonstrate parity rules

% Functions definition. h(x) should be odd, i(x) should be even
f = @(x) 2.*x.^2;
g = @(x) 3.*x.^3;
h = @(x) g(x).*f(x);
i = @(x) g(x).*g(x);

% Define the time domain. 400 divisions
x = linspace(-1, 1, 400);

% Simple mathematics [f(1)=2 g(1)=3 h(1)=6 i(x)=9]
disp([f(1), g(1), h(1), i(1)]);

% Plot
figure;
plot(x, f(x), 'LineWidth', 1.5); hold on;
plot(x, g(x), '--', 'LineWidth', 1.5);
plot(x, h(x),  ':', 'LineWidth', 1.5);
plot(x, i(x),  '.', 'LineWidth', 1.5);
grid on; xlabel('x'); ylabel('y');
legend('f(x) = 2x^2', 'g(x)=3x^3', 'h(x)=f(x).g(x)', 'i(x)=g(x).g(x)');
title('Plots of f, g, h, and i');