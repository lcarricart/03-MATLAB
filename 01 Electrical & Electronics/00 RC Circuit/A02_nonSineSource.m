% Purpose of the program: analysis of the compoment's voltages of an RC
% circuit using a <non-sinusoidal> source. 

% Problem: the phasors method assumes a single frequency, and an unchanging
% output when compared to the input signal. This works perfectly for
% sinusoidal inputs, but when using triangle waves or other very specific
% inputs, the phasors method is useless because the output signal doesn't
% resemble to the input.

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