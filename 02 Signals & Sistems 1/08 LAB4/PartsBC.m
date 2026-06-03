clc; clearvars;
% PART b
syms s t
Y = laplace((cos(2*t) + 2*sin(2*t))*heaviside(t), t, s);
pretty(Y)

% PART C
X = 1 + (s/(s^2 + 9));
x = ilaplace(X, s, t);
pretty(x)