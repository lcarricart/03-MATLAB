%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        %
%   SS2 Laboratory 2: Sampling Theorem - Problem 3       %
%                                                        %
%   Team Members: Georgii Molyboga  (2782258)            %
%                 Luciano Carricart (2782740)            %
%                 Mykyta  Kandyla   (2696614)            %
%                                                        %
%   Date: 03.06.2026                                     %
%                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the audio
[x_audio, fs_audio] = audioread('touchtone1.wav');

soundsc(x_audio, fs_audio); % Play sound

% Plot Time Domain
t_audio = (0:length(x_audio)-1) / fs_audio;

figure;
subplot(2,1,1);
plot(t_audio, x_audio);
xlabel('Time (s)');
title('Touchtone Signal in Time Domain');

% **********************************************************
%   STRATEGY 1 - Raw FFT
% **********************************************************
% Plot Frequency Domain
% X_audio = fft(x_audio);
% f_audio = (0:length(x_audio)-1) * (fs_audio/length(x_audio));

% **********************************************************
%   STRATEGY 2 - Actual Centered Spectrum
% **********************************************************
 N = length(x_audio);
 Ts = 1 / fs_audio;
 T_obs = N * Ts;
 Delta_f = 1 / T_obs;
 X_audio = fftshift(fft(x_audio)) * Ts;
 k = (-N/2 : N/2 - 1);          % (-floor(N/2) : ceil(N/2) - 1) could be even more precise
 f_audio = k * Delta_f;

subplot(2,1,2);
plot(f_audio, abs(X_audio));
xlabel('Frequency (Hz)');
title('Spectrum Over Entire Range');

% Zoomed in spectrum
figure;
plot(f_audio, abs(X_audio));
xlim([600 1600]);               % Range relevant for dual tones (~1455 max)
xlabel('Frequency (Hz)');
title('Zoomed Amplitude Spectrum');