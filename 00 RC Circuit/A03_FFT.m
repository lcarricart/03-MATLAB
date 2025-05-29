%% Purpose of the program: process the FFT of a set of data.

% In order to achieve this, I first need to clear the DC offset of the signal (clear the mean value).

%% Removes the mean value of the collected data for further analysis.
gpsAlt_mean = mean(gpsAltitude);
gpsAlt_noOffset = gpsAltitude - gpsAlt_mean;


% Clears Command Window & prints mean value
clc
fprintf("Mean value of gpsAltitude: %.2f \n", gpsAlt_mean)

% Plots the original signal without its offset
figure('Name','Original Plot','NumberTitle','off');
plot(time_ms, gpsAlt_noOffset, '.-', 'LineWidth', 0.6, 'MarkerSize', 1);
xlabel('Time [ms]');             
ylabel('GPS Altitude');
grid on;
title('Time series of sensor');

%% Fast Fourier Transform (Step 1: Uniform steps)

% We need a uniform sampling for the FFT (our current sampling goes around 1995 - 2010 points per row. We need this fixed).
delta_t  = 2000;                         % ms
t_start  = time_ms(1);                   % first timestamp
t_end    = time_ms(end);                 % last timestamp
t_uni    = t_start : delta_t : t_end;    % uniform grid in ms

gpsAlt_uni = interp1(      ...
    time_ms,               ...% original timestamps
    gpsAlt_noOffset,       ...% your altitude vector (DC removed)
    t_uni,                 ...% uniform timestamps
    'linear');             ...% interpolation method

% Plots the unified signal without its offset
figure('Name','Uniform Plot','NumberTitle','off');
plot(t_uni, gpsAlt_uni, '.-', 'LineWidth', 0.6, 'MarkerSize', 1, 'Color', "r");
xlabel('Time [ms] (uniform)');             
ylabel('GPS Altitude (uniform)');
grid on;
title('GPS Altitude: Uniform to 2000 ms / step');

% The signal looks correct to me :))

%% Fast Fourier Transform (Step 2: )

N   = numel(gpsAlt_uni);         % Number of points
fs  = 1 / 2000;                  % sampling frequency in Hz (fs = 0.5 for 2-s spacing)
% N2  = 2^nextpow2(N);           % optional padding
fourier   = fft(gpsAlt_uni, N);  % complex spectrum

P2  = abs(fourier).^2 / (fs*N);      % two-sided PSD (W/Hz or (unit)²/Hz)
P1  = P2(1:N/2);                        % one side PSD (keep positive-frequency half)
f   = fs*(0:N/2-1)/N;                   % 0 → fs/2  (0 → 0.25 Hz)

% Plots the FFT
figure
semilogx(f, P1)
grid on
xlabel('Frequency (Hz)')
ylabel('Power (W)')
title('One-sided FFT of uniformly-sampled GPS Altitude')