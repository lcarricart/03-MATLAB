% Propose of the program: process the full FFT anaylsis of one of the sensor's data.

% %% Data Import
% % READ Time
% time = readmatrix(                      ...
%                'Datalogger.xlsx',       ...
%                'Range','A:A',           ...      % first specific column
%                'OutputType','double');
% time = time(2:end);
% time = time - 2549585;                           % adjust the reference
% time_hours = time / 1000 / 60 / 60;
% time_ms = time;
% time_s = time_ms / 1000;
% 
% % READ GPS Altitude
% gpsAltitude = readmatrix(               ...
%                'Datalogger.xlsx',       ...
%                'Range','M:M',           ...      % second specific column
%                'OutputType','double');
% gpsAltitude = gpsAltitude(2:end);

%% Removes the mean value of the collected data for further analysis.
gpsAlt_mean = mean(gpsAltitude);
gpsAlt_noOffset = gpsAltitude - gpsAlt_mean;

% Clears Command Window & prints mean value
clc
fprintf("Mean value of gpsAltitude: %.2f \n", gpsAlt_mean)

% Plots the original signal without its offset
figure;
plot(time_s, gpsAlt_noOffset, '.-', 'LineWidth', 0.6, 'MarkerSize', 1);
xlabel('Time [s]');             
ylabel('GPS Altitude');
grid on;
title('Original Signal (no offset)');

%% Fast Fourier Transform (Step 1: Data pre-processing, uniform steps)

% We need a uniform sampling for the FFT (our current sampling goes around 1995 - 2010 points per row. We need this fixed).
delta_t  = 2;                            % seconds
t_start  = time_s(1);                    % first timestamp
t_end    = time_s(end);                  % last timestamp
t_uni    = t_start : delta_t : t_end;    % uniform grid in seconds

gpsAlt_uni = interp1(      ...
    time_s,                ...% original timestamps
    gpsAlt_noOffset,       ...% your altitude vector (DC removed)
    t_uni,                 ...% uniform timestamps
    'linear');             ...% interpolation method

% Plots the unified signal without its offset
figure;
plot(t_uni, gpsAlt_uni, '.-', 'LineWidth', 0.6, 'MarkerSize', 1, 'Color', "r");
xlabel('Time [s] (uniform)');             
ylabel('GPS Altitude (uniform)');
grid on;
title('Uniform signal, 2s / step');

% The signal looks correct to me :))

%% Fast Fourier Transform (Step 2: Perform the transform)

N   = numel(gpsAlt_uni);         % Number of points (lenght of the signal in terms of vector size)
Fs  = 1 / 2;                     % sampling frequency in Hz (fs = 0.5 for 2-s spacing)

% Compute the Fourier transform of the signal (Complex spectrum, real and imaginary)
fourier   = fft(gpsAlt_uni, N);

% Because Fourier transforms involve complex numbers, plot the complex magnitude of the fft spectrum.
figure;
plot(Fs/N*(0:N-1),abs(fourier),"LineWidth",3)
title("Complex Magnitude of fft Spectrum")
xlabel("f (Hz)")
ylabel("|fft(X)|")


%% Double-sided FFT Spectrum

%{ 
The plot shows five frequency peaks including the peak at 0 Hz for the DC offset. In this example, the signal is expected to have three frequency 
peaks at 0 Hz, 50 Hz, and 120 Hz. Here, the second half of the plot is the mirror reflection of the first half without including the peak at 0 Hz. 
The reason is that the discrete Fourier transform of a time-domain signal has a periodic nature, where the first half of its spectrum is in positive 
frequencies and the second half is in negative frequencies, with the first element reserved for the zero frequency.

For real signals, the fft spectrum is a two-sided spectrum, where the spectrum in the positive frequencies is the complex conjugate of the spectrum 
in the negative frequencies. To show the fft spectrum in the positive and negative frequencies, you can use fftshift. For an even length of L, the 
frequency domain starts from the negative of the Nyquist frequency -Fs/2 up to Fs/2-Fs/L with a spacing or frequency resolution of Fs/L.
%}

figure
plot(Fs/N*(-N/2:N/2-1),abs(fftshift(fourier)),"LineWidth",3)
title("fft Spectrum in the Positive and Negative Frequencies")
xlabel("f (Hz)")
ylabel("|fft(X)|")

%% Single-sided FFT Spectrum

%{
To find the amplitudes of the three frequency peaks, convert the fft spectrum in Y to the single-sided amplitude spectrum. Because the fft function 
includes a scaling factor L between the original and the transformed signals, rescale Y by dividing by L. Take the complex magnitude of the fft spectrum.
The two-sided amplitude spectrum P2, where the spectrum in the positive frequencies is the complex conjugate of the spectrum in the negative frequencies, 
has half the peak amplitudes of the time-domain signal. To convert to the single-sided spectrum, take the first half of the two-sided spectrum P2. Multiply 
the spectrum in the positive frequencies by 2. You do not need to multiply P1(1) and P1(end) by 2 because these amplitudes correspond to the zero and 
Nyquist frequencies, respectively, and they do not have the complex conjugate pairs in the negative frequencies.
%}

P2 = abs(fourier/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Define the frequency domain f for the single-sided spectrum. Plot the single-sided amplitude spectrum P1. As expected, the amplitudes are close to 0.8, 
% 0.7, and 1, but they are not exact because of the added noise. In most cases, longer signals produce better frequency approximations.
f = Fs/N*(0:(N/2));

figure
plot(f,P1,"LineWidth",3) 
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")

