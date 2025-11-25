%% Voltages in the series RC circuit

j = sqrt(-1);
V_amplitude = 5;        % V
R = 10;                 % Ohm
f = 500;            % frequency
T = 1 / f;              % Period

omega = 2*pi*f;         % rad/s
C = 1*10^(-3);          % F
Zc = 1/(j*omega*C);     % Ohm
t  = 0 : T/50 : 2*T;    % Time

% Circuit impedance
Z = R + Zc;

% Current
I = V_amplitude/Z;

% Voltage of the source
v_src = V_amplitude * cos(omega*t);

% Voltage across the resistor R
Vr = R*I;
v_r = abs(Vr) * cos(omega*t + angle(Vr));

% Voltage across the capacitor Zc
Vc = Zc*I;
v_c = abs(Vc) * cos(omega*t + angle(Vc));

subplot(1,1,1);
plot(t, v_src, t, v_r, t, v_c);
grid, title('Voltage measurement')
xlabel('t, s')
ylabel('V')
legend('V_{src}','V_r','V_c','Location', 'northeast')

% % Consider values V, Vr and Vc as the ends of the phasors. And to define
% % phasors enter their origins.
% 
% VV = [0 V];
% 
% VVc= [0 Vc];
% 
% VVr = [0 Vr];
% 
% % Plotting the voltage phasor V
% plot(real(VV), imag(VV))
% 
% % In order to view the phasors correcty the scale of the real axis should
% % be the same as that of the imaginary axis. This means a sqare plotting
% % frame which we obtain with:
% axis('square')
% 
% % Convinient axes produced by the statement
% axis([0 0.012 -0.006 0.006])
% 
% % We can use hold function to keep the graph and superimpose on in the
% % other two phasors, Vr and Vc
% hold on;
% plot(real(VVr),imag(VVr))
% plot(real(VVc),imag(VVc))
% 
% % Add the names
% title('Voltages in the series RC circuit')
% xlabel('Real'),ylabel('Imaginary')
% text(real(V), imag(V), 'V')
% text(real(Vr), imag(Vr), 'Vr')
% text(real(Vc), imag(Vc), 'Vc')
% 
% % To allow for future plots leave hold on state using hold off:
% hold off
% 
% %% Time plots of the claculated voltages
% 
% f = 50*10^6; % Frequency, Hz
% T = 1/f; % Period, s
% omega = 2*pi*f; % Angular frequency, rad/s
% t = 0:  T/50:  2*T; % Array of time values, s
% 
% v = V*sin(omega*t);
% vr = abs(Vr)*sin(omega*t + angle(Vr));
% vc = abs(Vc)*sin(omega*t + angle(Vc));
% 
% subplot(1, 2, 2);
% 
% plot(t, 1000*v, t, 1000*vr, t, 1000*vc)
% grid, title('Time plots of calculated voltages')
% xlabel('t,s'), ylabel('mV')
% text(t(5), 1000*v(5), 'V')
% text(t(20), 1000*v(20), 'Vr')
% text(t(50), 1000*v(50), 'Vc')
% 
% %% Magnitude and phase of the transfer function
% 
% fratio = 0:  0.01:  5;
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2,1,1), plot(fratio, abs(H))
% grid, ylabel('H(fratio)')
% 
% subplot(2,1,2), plot(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% 
% 
% %% Bode plot of the series RC lowpass filter
% 
% fratio = logspace(-1,2);
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2, 1, 1), semilogx(fratio, 20*log(abs(H)))
% grid, ylabel('Magnitude, dB')
% subplot(2, 1, 2), semilogx(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% % Consider values V, Vr and Vc as the ends of the phasors. And to define
% % phasors enter their origins.
% 
% VV = [0 V];
% 
% VVc= [0 Vc];
% 
% VVr = [0 Vr];
% 
% % Plotting the voltage phasor V
% plot(real(VV), imag(VV))
% 
% % In order to view the phasors correcty the scale of the real axis should
% % be the same as that of the imaginary axis. This means a sqare plotting
% % frame which we obtain with:
% axis('square')
% 
% % Convinient axes produced by the statement
% axis([0 0.012 -0.006 0.006])
% 
% % We can use hold function to keep the graph and superimpose on in the
% % other two phasors, Vr and Vc
% hold on;
% plot(real(VVr),imag(VVr))
% plot(real(VVc),imag(VVc))
% 
% % Add the names
% title('Voltages in the series RC circuit')
% xlabel('Real'),ylabel('Imaginary')
% text(real(V), imag(V), 'V')
% text(real(Vr), imag(Vr), 'Vr')
% text(real(Vc), imag(Vc), 'Vc')
% 
% % To allow for future plots leave hold on state using hold off:
% hold off
% 
% %% Time plots of the claculated voltages
% 
% f = 50*10^6; % Frequency, Hz
% T = 1/f; % Period, s
% omega = 2*pi*f; % Angular frequency, rad/s
% t = 0:  T/50:  2*T; % Array of time values, s
% 
% v = V*sin(omega*t);
% vr = abs(Vr)*sin(omega*t + angle(Vr));
% vc = abs(Vc)*sin(omega*t + angle(Vc));
% 
% subplot(1, 2, 2);
% 
% plot(t, 1000*v, t, 1000*vr, t, 1000*vc)
% grid, title('Time plots of calculated voltages')
% xlabel('t,s'), ylabel('mV')
% text(t(5), 1000*v(5), 'V')
% text(t(20), 1000*v(20), 'Vr')
% text(t(50), 1000*v(50), 'Vc')
% 
% %% Magnitude and phase of the transfer function
% 
% fratio = 0:  0.01:  5;
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2,1,1), plot(fratio, abs(H))
% grid, ylabel('H(fratio)')
% 
% subplot(2,1,2), plot(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% 
% 
% %% Bode plot of the series RC lowpass filter
% 
% fratio = logspace(-1,2);
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2, 1, 1), semilogx(fratio, 20*log(abs(H)))
% grid, ylabel('Magnitude, dB')
% subplot(2, 1, 2), semilogx(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% % Consider values V, Vr and Vc as the ends of the phasors. And to define
% % phasors enter their origins.
% 
% VV = [0 V];
% 
% VVc= [0 Vc];
% 
% VVr = [0 Vr];
% 
% % Plotting the voltage phasor V
% plot(real(VV), imag(VV))
% 
% % In order to view the phasors correcty the scale of the real axis should
% % be the same as that of the imaginary axis. This means a sqare plotting
% % frame which we obtain with:
% axis('square')
% 
% % Convinient axes produced by the statement
% axis([0 0.012 -0.006 0.006])
% 
% % We can use hold function to keep the graph and superimpose on in the
% % other two phasors, Vr and Vc
% hold on;
% plot(real(VVr),imag(VVr))
% plot(real(VVc),imag(VVc))
% 
% % Add the names
% title('Voltages in the series RC circuit')
% xlabel('Real'),ylabel('Imaginary')
% text(real(V), imag(V), 'V')
% text(real(Vr), imag(Vr), 'Vr')
% text(real(Vc), imag(Vc), 'Vc')
% 
% % To allow for future plots leave hold on state using hold off:
% hold off
% 
% %% Time plots of the claculated voltages
% 
% f = 50*10^6; % Frequency, Hz
% T = 1/f; % Period, s
% omega = 2*pi*f; % Angular frequency, rad/s
% t = 0:  T/50:  2*T; % Array of time values, s
% 
% v = V*sin(omega*t);
% vr = abs(Vr)*sin(omega*t + angle(Vr));
% vc = abs(Vc)*sin(omega*t + angle(Vc));
% 
% subplot(1, 2, 2);
% 
% plot(t, 1000*v, t, 1000*vr, t, 1000*vc)
% grid, title('Time plots of calculated voltages')
% xlabel('t,s'), ylabel('mV')
% text(t(5), 1000*v(5), 'V')
% text(t(20), 1000*v(20), 'Vr')
% text(t(50), 1000*v(50), 'Vc')
% 
% %% Magnitude and phase of the transfer function
% 
% fratio = 0:  0.01:  5;
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2,1,1), plot(fratio, abs(H))
% grid, ylabel('H(fratio)')
% 
% subplot(2,1,2), plot(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% 
% 
% %% Bode plot of the series RC lowpass filter
% 
% fratio = logspace(-1,2);
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2, 1, 1), semilogx(fratio, 20*log(abs(H)))
% grid, ylabel('Magnitude, dB')
% subplot(2, 1, 2), semilogx(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% % Consider values V, Vr and Vc as the ends of the phasors. And to define
% % phasors enter their origins.
% 
% VV = [0 V];
% 
% VVc= [0 Vc];
% 
% VVr = [0 Vr];
% 
% % Plotting the voltage phasor V
% plot(real(VV), imag(VV))
% 
% % In order to view the phasors correcty the scale of the real axis should
% % be the same as that of the imaginary axis. This means a sqare plotting
% % frame which we obtain with:
% axis('square')
% 
% % Convinient axes produced by the statement
% axis([0 0.012 -0.006 0.006])
% 
% % We can use hold function to keep the graph and superimpose on in the
% % other two phasors, Vr and Vc
% hold on;
% plot(real(VVr),imag(VVr))
% plot(real(VVc),imag(VVc))
% 
% % Add the names
% title('Voltages in the series RC circuit')
% xlabel('Real'),ylabel('Imaginary')
% text(real(V), imag(V), 'V')
% text(real(Vr), imag(Vr), 'Vr')
% text(real(Vc), imag(Vc), 'Vc')
% 
% % To allow for future plots leave hold on state using hold off:
% hold off
% 
% %% Time plots of the claculated voltages
% 
% f = 50*10^6; % Frequency, Hz
% T = 1/f; % Period, s
% omega = 2*pi*f; % Angular frequency, rad/s
% t = 0:  T/50:  2*T; % Array of time values, s
% 
% v = V*sin(omega*t);
% vr = abs(Vr)*sin(omega*t + angle(Vr));
% vc = abs(Vc)*sin(omega*t + angle(Vc));
% 
% subplot(1, 2, 2);
% 
% plot(t, 1000*v, t, 1000*vr, t, 1000*vc)
% grid, title('Time plots of calculated voltages')
% xlabel('t,s'), ylabel('mV')
% text(t(5), 1000*v(5), 'V')
% text(t(20), 1000*v(20), 'Vr')
% text(t(50), 1000*v(50), 'Vc')
% 
% %% Magnitude and phase of the transfer function
% 
% fratio = 0:  0.01:  5;
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2,1,1), plot(fratio, abs(H))
% grid, ylabel('H(fratio)')
% 
% subplot(2,1,2), plot(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% 
% 
% %% Bode plot of the series RC lowpass filter
% 
% fratio = logspace(-1,2);
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2, 1, 1), semilogx(fratio, 20*log(abs(H)))
% grid, ylabel('Magnitude, dB')
% subplot(2, 1, 2), semilogx(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% % Consider values V, Vr and Vc as the ends of the phasors. And to define
% % phasors enter their origins.
% 
% VV = [0 V];
% 
% VVc= [0 Vc];
% 
% VVr = [0 Vr];
% 
% % Plotting the voltage phasor V
% plot(real(VV), imag(VV))
% 
% % In order to view the phasors correcty the scale of the real axis should
% % be the same as that of the imaginary axis. This means a sqare plotting
% % frame which we obtain with:
% axis('square')
% 
% % Convinient axes produced by the statement
% axis([0 0.012 -0.006 0.006])
% 
% % We can use hold function to keep the graph and superimpose on in the
% % other two phasors, Vr and Vc
% hold on;
% plot(real(VVr),imag(VVr))
% plot(real(VVc),imag(VVc))
% 
% % Add the names
% title('Voltages in the series RC circuit')
% xlabel('Real'),ylabel('Imaginary')
% text(real(V), imag(V), 'V')
% text(real(Vr), imag(Vr), 'Vr')
% text(real(Vc), imag(Vc), 'Vc')
% 
% % To allow for future plots leave hold on state using hold off:
% hold off
% 
% %% Time plots of the claculated voltages
% 
% f = 50*10^6; % Frequency, Hz
% T = 1/f; % Period, s
% omega = 2*pi*f; % Angular frequency, rad/s
% t = 0:  T/50:  2*T; % Array of time values, s
% 
% v = V*sin(omega*t);
% vr = abs(Vr)*sin(omega*t + angle(Vr));
% vc = abs(Vc)*sin(omega*t + angle(Vc));
% 
% subplot(1, 2, 2);
% 
% plot(t, 1000*v, t, 1000*vr, t, 1000*vc)
% grid, title('Time plots of calculated voltages')
% xlabel('t,s'), ylabel('mV')
% text(t(5), 1000*v(5), 'V')
% text(t(20), 1000*v(20), 'Vr')
% text(t(50), 1000*v(50), 'Vc')
% 
% %% Magnitude and phase of the transfer function
% 
% fratio = 0:  0.01:  5;
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2,1,1), plot(fratio, abs(H))
% grid, ylabel('H(fratio)')
% 
% subplot(2,1,2), plot(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% 
% 
% %% Bode plot of the series RC lowpass filter
% 
% fratio = logspace(-1,2);
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2, 1, 1), semilogx(fratio, 20*log(abs(H)))
% grid, ylabel('Magnitude, dB')
% subplot(2, 1, 2), semilogx(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% % Consider values V, Vr and Vc as the ends of the phasors. And to define
% % phasors enter their origins.
% 
% VV = [0 V];
% 
% VVc= [0 Vc];
% 
% VVr = [0 Vr];
% 
% % Plotting the voltage phasor V
% plot(real(VV), imag(VV))
% 
% % In order to view the phasors correcty the scale of the real axis should
% % be the same as that of the imaginary axis. This means a sqare plotting
% % frame which we obtain with:
% axis('square')
% 
% % Convinient axes produced by the statement
% axis([0 0.012 -0.006 0.006])
% 
% % We can use hold function to keep the graph and superimpose on in the
% % other two phasors, Vr and Vc
% hold on;
% plot(real(VVr),imag(VVr))
% plot(real(VVc),imag(VVc))
% 
% % Add the names
% title('Voltages in the series RC circuit')
% xlabel('Real'),ylabel('Imaginary')
% text(real(V), imag(V), 'V')
% text(real(Vr), imag(Vr), 'Vr')
% text(real(Vc), imag(Vc), 'Vc')
% 
% % To allow for future plots leave hold on state using hold off:
% hold off
% 
% %% Time plots of the claculated voltages
% 
% f = 50*10^6; % Frequency, Hz
% T = 1/f; % Period, s
% omega = 2*pi*f; % Angular frequency, rad/s
% t = 0:  T/50:  2*T; % Array of time values, s
% 
% v = V*sin(omega*t);
% vr = abs(Vr)*sin(omega*t + angle(Vr));
% vc = abs(Vc)*sin(omega*t + angle(Vc));
% 
% subplot(1, 2, 2);
% 
% plot(t, 1000*v, t, 1000*vr, t, 1000*vc)
% grid, title('Time plots of calculated voltages')
% xlabel('t,s'), ylabel('mV')
% text(t(5), 1000*v(5), 'V')
% text(t(20), 1000*v(20), 'Vr')
% text(t(50), 1000*v(50), 'Vc')
% 
% %% Magnitude and phase of the transfer function
% 
% fratio = 0:  0.01:  5;
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2,1,1), plot(fratio, abs(H))
% grid, ylabel('H(fratio)')
% 
% subplot(2,1,2), plot(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% 
% 
% %% Bode plot of the series RC lowpass filter
% 
% fratio = logspace(-1,2);
% H = ones(size(fratio))./(1 + j*fratio);
% 
% subplot(2, 1, 1), semilogx(fratio, 20*log(abs(H)))
% grid, ylabel('Magnitude, dB')
% subplot(2, 1, 2), semilogx(fratio, 180*angle(H)/pi)
% grid, xlabel('f/f_c'), ylabel('Phase, deg')
% Consider values V, Vr and Vc as the ends of the phasors. And to define
% phasors enter their origins.

VV = [0 V];

VVc= [0 Vc];

VVr = [0 Vr];

% Plotting the voltage phasor V
plot(real(VV), imag(VV))

% In order to view the phasors correcty the scale of the real axis should
% be the same as that of the imaginary axis. This means a sqare plotting
% frame which we obtain with:
axis('square')

% Convinient axes produced by the statement
axis([0 0.012 -0.006 0.006])

% We can use hold function to keep the graph and superimpose on in the
% other two phasors, Vr and Vc
hold on;
plot(real(VVr),imag(VVr))
plot(real(VVc),imag(VVc))

% Add the names
title('Voltages in the series RC circuit')
xlabel('Real'),ylabel('Imaginary')
text(real(V), imag(V), 'V')
text(real(Vr), imag(Vr), 'Vr')
text(real(Vc), imag(Vc), 'Vc')

% To allow for future plots leave hold on state using hold off:
hold off

%% Time plots of the claculated voltages

f = 50*10^6; % Frequency, Hz
T = 1/f; % Period, s
omega = 2*pi*f; % Angular frequency, rad/s
t = 0:  T/50:  2*T; % Array of time values, s

v = V*sin(omega*t);
vr = abs(Vr)*sin(omega*t + angle(Vr));
vc = abs(Vc)*sin(omega*t + angle(Vc));

subplot(1, 2, 2);

plot(t, 1000*v, t, 1000*vr, t, 1000*vc)
grid, title('Time plots of calculated voltages')
xlabel('t,s'), ylabel('mV')
text(t(5), 1000*v(5), 'V')
text(t(20), 1000*v(20), 'Vr')
text(t(50), 1000*v(50), 'Vc')

%% Magnitude and phase of the transfer function

fratio = 0:  0.01:  5;
H = ones(size(fratio))./(1 + j*fratio);

subplot(2,1,1), plot(fratio, abs(H))
grid, ylabel('H(fratio)')

subplot(2,1,2), plot(fratio, 180*angle(H)/pi)
grid, xlabel('f/f_c'), ylabel('Phase, deg')


%% Bode plot of the series RC lowpass filter

fratio = logspace(-1,2);
H = ones(size(fratio))./(1 + j*fratio);

subplot(2, 1, 1), semilogx(fratio, 20*log(abs(H)))
grid, ylabel('Magnitude, dB')
subplot(2, 1, 2), semilogx(fratio, 180*angle(H)/pi)
grid, xlabel('f/f_c'), ylabel('Phase, deg')
% Consider values V, Vr and Vc as the ends of the phasors. And to define
% phasors enter their origins.

VV = [0 V];

VVc= [0 Vc];

VVr = [0 Vr];

% Plotting the voltage phasor V
plot(real(VV), imag(VV))

% In order to view the phasors correcty the scale of the real axis should
% be the same as that of the imaginary axis. This means a sqare plotting
% frame which we obtain with:
axis('square')

% Convinient axes produced by the statement
axis([0 0.012 -0.006 0.006])

% We can use hold function to keep the graph and superimpose on in the
% other two phasors, Vr and Vc
hold on;
plot(real(VVr),imag(VVr))
plot(real(VVc),imag(VVc))

% Add the names
title('Voltages in the series RC circuit')
xlabel('Real'),ylabel('Imaginary')
text(real(V), imag(V), 'V')
text(real(Vr), imag(Vr), 'Vr')
text(real(Vc), imag(Vc), 'Vc')

% To allow for future plots leave hold on state using hold off:
hold off

%% Time plots of the claculated voltages

f = 50*10^6; % Frequency, Hz
T = 1/f; % Period, s
omega = 2*pi*f; % Angular frequency, rad/s
t = 0:  T/50:  2*T; % Array of time values, s

v = V*sin(omega*t);
vr = abs(Vr)*sin(omega*t + angle(Vr));
vc = abs(Vc)*sin(omega*t + angle(Vc));

subplot(1, 2, 2);

plot(t, 1000*v, t, 1000*vr, t, 1000*vc)
grid, title('Time plots of calculated voltages')
xlabel('t,s'), ylabel('mV')
text(t(5), 1000*v(5), 'V')
text(t(20), 1000*v(20), 'Vr')
text(t(50), 1000*v(50), 'Vc')

%% Magnitude and phase of the transfer function

fratio = 0:  0.01:  5;
H = ones(size(fratio))./(1 + j*fratio);

subplot(2,1,1), plot(fratio, abs(H))
grid, ylabel('H(fratio)')

subplot(2,1,2), plot(fratio, 180*angle(H)/pi)
grid, xlabel('f/f_c'), ylabel('Phase, deg')


%% Bode plot of the series RC lowpass filter

fratio = logspace(-1,2);
H = ones(size(fratio))./(1 + j*fratio);

subplot(2, 1, 1), semilogx(fratio, 20*log(abs(H)))
grid, ylabel('Magnitude, dB')
subplot(2, 1, 2), semilogx(fratio, 180*angle(H)/pi)
grid, xlabel('f/f_c'), ylabel('Phase, deg')
% Consider values V, Vr and Vc as the ends of the phasors. And to define
% phasors enter their origins.

VV = [0 V];

VVc= [0 Vc];

VVr = [0 Vr];

% Plotting the voltage phasor V
plot(real(VV), imag(VV))

% In order to view the phasors correcty the scale of the real axis should
% be the same as that of the imaginary axis. This means a sqare plotting
% frame which we obtain with:
axis('square')

% Convinient axes produced by the statement
axis([0 0.012 -0.006 0.006])

% We can use hold function to keep the graph and superimpose on in the
% other two phasors, Vr and Vc
hold on;
plot(real(VVr),imag(VVr))
plot(real(VVc),imag(VVc))

% Add the names
title('Voltages in the series RC circuit')
xlabel('Real'),ylabel('Imaginary')
text(real(V), imag(V), 'V')
text(real(Vr), imag(Vr), 'Vr')
text(real(Vc), imag(Vc), 'Vc')

% To allow for future plots leave hold on state using hold off:
hold off

%% Time plots of the claculated voltages

f = 50*10^6; % Frequency, Hz
T = 1/f; % Period, s
omega = 2*pi*f; % Angular frequency, rad/s
t = 0:  T/50:  2*T; % Array of time values, s

v = V*sin(omega*t);
vr = abs(Vr)*sin(omega*t + angle(Vr));
vc = abs(Vc)*sin(omega*t + angle(Vc));

subplot(1, 2, 2);

plot(t, 1000*v, t, 1000*vr, t, 1000*vc)
grid, title('Time plots of calculated voltages')
xlabel('t,s'), ylabel('mV')
text(t(5), 1000*v(5), 'V')
text(t(20), 1000*v(20), 'Vr')
text(t(50), 1000*v(50), 'Vc')

%% Magnitude and phase of the transfer function

fratio = 0:  0.01:  5;
H = ones(size(fratio))./(1 + j*fratio);

subplot(2,1,1), plot(fratio, abs(H))
grid, ylabel('H(fratio)')

subplot(2,1,2), plot(fratio, 180*angle(H)/pi)
grid, xlabel('f/f_c'), ylabel('Phase, deg')


%% Bode plot of the series RC lowpass filter

fratio = logspace(-1,2);
H = ones(size(fratio))./(1 + j*fratio);

subplot(2, 1, 1), semilogx(fratio, 20*log(abs(H)))
grid, ylabel('Magnitude, dB')
subplot(2, 1, 2), semilogx(fratio, 180*angle(H)/pi)
grid, xlabel('f/f_c'), ylabel('Phase, deg')
