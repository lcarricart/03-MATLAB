function [a_k, b_k, x_fourier] = fourier_square_wave(N, f0, t)
% FOURIER_SQUARE_WAVE Computes Fourier series of a square wave
%   Inputs:
%       N  - Number of Fourier series coefficients to compute
%       f0 - Fundamental frequency of the square wave [Hz]
%       t  - Time vector for computing the Fourier series
%
%   Outputs:
%       a_k, b_k fourier coefficients
%       x_fourier - Reconstructed signal using Fourier series

    % Initialize coefficients
    a_k = zeros(1, N+1);  % a_0, a_1, ..., a_N
    b_k = zeros(1, N);    % b_1, ..., b_N
    
    % For square wave toggling between -1 and +1:
    % a_0 = 0 (DC component is zero)
    % a_k = 0 for all k (even function property)
    % b_k = 4/(pi*k) for odd k, 0 for even k
    
    a_k(1) = 0;  % a_0 = 0
    
    for k = 1:N
        a_k(k+1) = 0;  % All cosine coefficients are 0
        
        if mod(k, 2) == 1  % Odd harmonics only
            b_k(k) = 4 / (pi * k);
        else
            b_k(k) = 0;
        end
    end
    
    % Reconstruct signal using Fourier series
    % x(t) = a_0 + sum(a_k*cos(2*pi*k*f0*t) + b_k*sin(2*pi*k*f0*t))
    x_fourier = a_k(1) * ones(size(t));  % DC component
    
    for k = 1:N
        x_fourier = x_fourier + a_k(k+1) * cos(2*pi*k*f0*t) + b_k(k) * sin(2*pi*k*f0*t);
    end
    
end
