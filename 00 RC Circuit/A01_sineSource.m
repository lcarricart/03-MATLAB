% Purpose of the program: analysis of the compoment's voltages of an RC
% circuit using a sinusoidal source.

%% Voltages in the series RC circuit

j = sqrt(-1);
V_amplitude = 5;        % V
R = 10;                 % Ohm
f = 500;            % frequency
T = 1 / f;              % Period

omega = 2*pi*f;         % rad/s
C = 1*10^(-3);          % F
Zc = 1/(j*omega*C);     % Ohm
t  = 0 : T/50 : 2*T;    % Time (from 0 to 2 periods, with a step of T/50

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