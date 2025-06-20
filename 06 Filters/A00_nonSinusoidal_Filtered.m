% RCFilter
clear all;

R = 100; C = 4.7e-6;           % Resistor, Capacitor
n = 49;                        % number of harmonics
f0 = 100;                      % Frequency in Hz
t  = [0:0.01:20]*1e-3;         % Time in s

for k = 1:2:n
    % Fourier coefficients of square signal
    % A(k)=4/(pi*k) if k is uneven, else 0; phase always 0
    uAmpIn(k)  = 4/(pi*k);
    uin(k,:)   = uAmpIn(k).*sin(2*pi*k*f0.*t);

    % Transfer function at k*2*pi*f
    F          = 1/(1 + j*2*pi*k*f0*R*C);
    uAmpOut(k) = uAmpIn(k)*F;
    uOut(k,:)  = imag(uAmpOut(k).*exp(j*2*pi*k*f0.*t));
end

figure(1);
clf;
plot(t*1000, sum(uin),  'g-', 'LineWidth', 2);   % Time axis in ms
hold on;
plot(t*1000, sum(uOut), 'r-', 'LineWidth', 2);   % Time axis in ms
hold off;

ylim([-2 2]);
xlabel('t / ms');
ylabel('u / V');
title('RC lowpass filtering of sinusoidal signal');
text(17, 1.7, ['n = ', num2str(k)]);
grid on;
legend('input', 'output');
legend('Location', 'southeast');
