%% Purpose of the program: compute the response of an RLC circuit. Series RLC in parallel with L. 
% L2.i''' + L1.i'' + R.i' + C1.i = 0

% Circuit's parameters
R  = 100;         % Ohm
L1 = 1e-3;     % H 
L2 = 100e-3;      % H
C  = 47e-3;     % F

% Define the ODE function
odefun = @(t, y) [ y(2);
                   y(3);
                   ( -L1*y(3) - R*y(2) - (1/C)*y(1) ) / L2 ];

% Define the time span and the initial conditions [i(0); i'(0); i''(0)]
tspan = [0, 100];
y0 = [10; 0; -100];

% Call MATLAB's ODE solver
[t, y] = ode45(odefun, tspan, y0);

% Plot
plot(t, y(:,1), 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Current i(t) (A)')
title('Third-order Series RLC Circuit with Two Inductors')
grid on


% This answer represents an unstable circuit. This is unreproduceable because of the initial condition i''(0) = -100 A/s^2, so it is in fact a
% bit unreal. However, this shows that for this value of R, the circuit has a positive pole.