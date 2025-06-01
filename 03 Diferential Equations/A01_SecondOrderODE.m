%% Purpose of the program: solve the second order ODE y'' + 2y' + y = 0 with the initial conditions y(0) = 1, and y'(0) = 0

% Define the ODE function:
%   - Quick way for simple functions (used in my program). @(t, y)
%   - Complex functions: define a separate function file describing the equation

% MATLAB doesn't see DE the same way humans do, so if your ODE is second order (involves y''), MATLAB needs you to break it into:
% “How does y  change?” → (y')
% “How does y' change?” → (y'' written in terms of y and y')
% y(1) is y
% y(2) is y'

odefun = @(t, y) [y(2);
                 -2*y(2) - y(1)];

% Define the time span and the initial conditions
tspan = [0 10];              % Time from 0 to 10
y0 = [1; 0];                 % Initial conditions: y(0)=1, y'(0)=0

% Call MATLAB's ODE solver
[t, y] = ode45(odefun, tspan, y0);

plot(t, y(:,1));             % y(:,1) is y(t)
xlabel('Time t');
ylabel('y(t)');
title('Solution of y'''' + 2y'' + y = 0');
grid on;
