function [x_e, x_o] = decompose_signal(x)
% DECOMPOSE_SIGNAL Decomposes a signal into even and odd parts
%   Input:
%       x - vector of length N (odd number)
%           Format: [(N-1)/2 values for t<0], [value at t=0], [(N-1)/2 values for t>0]
%   
%   Output:
%       x_e - even part of x
%       x_o - odd part of x
%
%   Formulas:
%       x_e(t) = 0.5 * [x(t) + x(-t)]
%       x_o(t) = 0.5 * [x(t) - x(-t)]

    % Get length of signal
    N = length(x);
    
    % Create time-reversed version of x (flip the signal)
    x_rev = fliplr(x);
    
    % Calculate even part: x_e(t) = 0.5 * [x(t) + x(-t)]
    x_e = 0.5 * (x + x_rev);
    
    % Calculate odd part: x_o(t) = 0.5 * [x(t) - x(-t)]
    x_o = 0.5 * (x - x_rev);
    
end
