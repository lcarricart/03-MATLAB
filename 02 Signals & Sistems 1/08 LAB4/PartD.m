clc; clearvars;

syms t s y(t) Y

% Define ODE
ode = diff(y, t, 2) + 3*diff(y, t) + 2*y == exp(-t)*heaviside(t);

% Laplace transform
L_ode = laplace(ode, t, s)

% Substitute Laplace{y(t)} = Y
L_ode = subs(L_ode, laplace(y(t), t, s), Y);

% Apply initial conditions
L_ode = subs(L_ode, {y(0), subs(diff(y,t),t,0)}, {0,0});

% Solve for Y(s)
Y = solve(L_ode, Y);
Y = simplify(Y)

% Inverse Laplace
y = ilaplace(Y, s, t);
pretty(y)