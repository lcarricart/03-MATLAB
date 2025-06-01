%% Purpose of the program: analyze the oscillating RLC circuit proposed for Electrical Engineering 2, L.y'' + R.y' + (1/C).y = 0 

% Circuit's parameters
R = 1;
L = 10e-3;
C = 47e-6;

% Define the ODE function
odefun = @(t, y) [y(2); (-R/L)*y(2) - (1/(L*C))*y(1)];

% Define the time span and the initial conditions
tspan = [0 0.1];               
y0 = [0; -100];                 % Initial conditions: y(0)=0, y'(0)= -100

% Call MATLAB's ODE solver
[t, y] = ode45(odefun, tspan, y0);

plot(t, y(:,1));             % y(:,1) is y(t)
xlabel('Time t');
ylabel('y(t)');
title('Solution');
grid on;
