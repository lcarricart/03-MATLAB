%% Purpose of the program: analyze how a FFT works with a simple function.

%% Parameters for a clean test
fs  = 1000;          % sampling rate, Hz
f0  = 123;           % sine-wave frequency, Hz  (any value < fs/2)
N   = 1024;          % number of samples (power of two for convenience)

%% 1  Build a uniformly-sampled sine wave
t = (0:N-1)/fs;      % bulding the time vector, seconds
x = sin(2*pi*f0*t);  % 1-k amplitude sine

%% 2  Compute the FFT (no window, no zero-padding)
X = fft(x);          % complex spectrum, N points

% %% 3  Extract the one-sided amplitude spectrum
% % Fourier Transforms of a sine wave would show two peaks, of equal value but opposite signs. It`s mostly insteresting for us to see one of them
% P = abs(X)/N;        % scale by N so peak = 0.5
% P = 2*P(1:N/2);      % keeps half of the spectrum (N/2) and restore full amplitude that was lost in the last step
% f = fs*(0:N/2-1)/N;  % matching frequency axis, 0 → Nyquist

%% 4  Plot
figure
stem(f, P, 'filled')
grid on
xlabel('frequency (Hz)')
ylabel('amplitude')
title(sprintf('one-sided FFT of %.0f Hz sine', f0))
xlim([0 fs/2])
