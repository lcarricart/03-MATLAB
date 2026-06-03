%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        %
%   SS2 Laboratory 2: Sampling Theorem - Problem 2       %
%                                                        %
%   Team Members: Georgii Molyboga  (2782258)            %
%                 Luciano Carricart (2782740)            %
%                 Mykyta  Kandyla   (2696614)            %
%                                                        %
%   Date: 03.06.2026                                     %
%                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 8000;          % Sampling frequency 8kHz
Ts = 1 / fs;
N = 64;             % Power of 2, 64 samples at 8kHz = 8ms  

t = (0:N-1) * Ts;
x = zeros(1, N);
x(t <= 0.002) = 1;  % Rectangular pulse

subplot(2,1,1);
stem(t*1000, x);
xlabel('Time (ms)');
ylabel('Amplitude');
title('Sampled Rectangular Signal');

% **********************************************************
%   STRATEGY 1 - Custom myDFT()
% **********************************************************
% Calculate spectrum using the custom myDFT function
% X = myDFT(x);
% T_obs = N * Ts;
% Delta_f = 1 / T_obs;              % Frequency resolution, also called frequency-bin spacing
% k = (-(N-1)/2 : (N-1)/2);
% f = k * Delta_f;                  % k is the bin index. For k=0, f=0. For k=1, f=Delta_f, ...

% **********************************************************
%   STRATEGY 2 - MATLAB's fft()
% **********************************************************
 X = fftshift(fft(x)) * Ts;
 T_obs = N * Ts;
 Delta_f = 1 / T_obs;
 k = (-N/2 : N/2 - 1);              % k is adjusted to use an even N with fftshift
 f = k * Delta_f;

subplot(2,1,2);
positive_idx = (f >= 0);            % Filter for positive frequencies to plot the range [0, fs/2]

% Plot the amplitude spectrum 
plot(f(positive_idx), abs(X(positive_idx)));
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
title('Amplitude Spectrum');