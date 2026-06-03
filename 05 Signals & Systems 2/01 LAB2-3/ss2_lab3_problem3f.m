%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        %
%   SS2 Laboratory 2: Sampling Theorem - Problem 3(f)    %
%                                                        %
%   Team Members: Georgii Molyboga  (2782258)            %
%                 Luciano Carricart (2782740)            %
%                 Mykyta  Kandyla   (2696614)            %
%                                                        %
%   Date: 07.05.2026                                     %
%                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 8000;                  % Sampling frequency (8 kHz)
my_number = '017612345678';

y = generateDTMF(my_number);
soundsc(y, fs);

% Time Domain Plot
t = (0:length(y)-1) / fs;

subplot(2,1,1);
plot(t, y);
xlabel('Time (s)');
ylabel('Amplitude');
title(['Time Domain: DTMF Signal for ', my_number]);

% Frequency Domain Plot
Y_fft = fft(y);
N = length(y);
f = (0:N-1) * (fs / N);

% We only need the positive frequencies up to the Nyquist limit (fs/2)
half_idx = floor(N/2) + 1;
f_pos = f(1:half_idx);
Y_mag = abs(Y_fft(1:half_idx));

subplot(2,1,2);
plot(f_pos, Y_mag);
% Restrict the x-axis to the relevant DTMF frequency band (600 Hz - 1600 Hz)
xlim([600, 1600]); 
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Domain: Amplitude Spectrum (Zoomed to DTMF Band)');
grid on;