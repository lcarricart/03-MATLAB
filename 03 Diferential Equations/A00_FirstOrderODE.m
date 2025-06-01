%% Purpose of the program: solve the first order ODE y' = -2y, with y(0) = 1

% Define the ODE as a function. For simple cases like first order ODEs, you can use the matlab anonymous function. Otherwise, create a separate
% function file.
odefun = @(t, y) -2*y;

% Define the time span and initial condition
tspan = [0 5];          
y0 = 1;                  

% Call the MATLAB's ODE solver
[t, y] = ode45(odefun, tspan, y0);

plot(t, y, 'LineWidth', 2)
xlabel('Time t')
ylabel('y(t)')
title('Solution of dy/dt = -2y with y(0) = 1')
grid on

% The solution shows an exponentially decaying function, starting at 1 (because of my initial condition y(0) = 1). This graph represents the
% behavior of the system and indicates that the solution has the form of y(0) * e^(-a t)