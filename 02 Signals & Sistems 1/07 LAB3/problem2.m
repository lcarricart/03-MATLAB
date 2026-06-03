%% Start of the script
clc;
define_plot_X_H();
define_plot_X_H2();
define_plot_X_H3();

%% Functions

% Discrete convolution of the original function x(t) and an impulse response h(t). The function can be used generally to convolve two causal functions.
% Input: x (original function OR function 1)
% Output: h (impulse response OR function 2)
function y = myConv(x, h)
    % In order to perform a convolution, it is not necessary that both vectors are of equal lenght.
    % A discrete convolution defines "y" for specific values of "t"
    lengthX = length(x);
    lengthH = length(h);
    
    % n will be the number of "y" values that I will get. It is the sum of length1 + length 2 because the possible values span from the first touch of the functions until the last touch of the functions. All other values are 0.
    n = lengthH + lengthX - 1;
    y = zeros(1, n);                        % Initialize the output vector y with zeros
   
    for n2 = 1:n
        sum = 0;
        
        for k = 1:lengthX
            m = n2 - k + 1;                  % This corresponds to index of h[n-k]
            
            if m >= 1 && m <= lengthH       % 1) m cant be smaller than 1, because a MATLAB array starts indexing from 1. 2) m can't be greater than lengthH because the indexing h(m) would cause an error.
                sum = sum + x(k) * h(m);
            end
        end

        y(n2) = sum;  % Store the computed sum in the output vector "y", position n2.
    end
end

function define_plot_X_H()
    x = [0, 1, 1, 1, 1, 0, 0, 0, 0, 0];
    h = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0];

    y = myConv(x, h);

    figure(Name='Exercise 2)B)');
    subplot(2, 1, 1);
    stem(x, 'filled', 'LineStyle', 'none', 'MarkerSize', 8);
    title('Input x(t) of the Convolution');
    xlabel('Time (s)');
    ylabel('Magnitude');
    xline(0);             % vertical line at ω = 0
    yline(0);             % horizontal line at 0 magnitude

    subplot(2, 1, 2);
    stem(y, 'filled', 'LineStyle', 'none', 'MarkerSize', 8);
    title('Output y(t) of the Convolution for h = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]');
    xlabel('Time (s)');
    ylabel('Magnitude');
    xline(0);             
    yline(0);             
end

function define_plot_X_H2()
    x = [0, 1, 1, 1, 1, 0, 0, 0, 0, 0];
    h = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0];

    y = myConv(x, h);

    figure(Name='Exercise 2)C)');
    subplot(2, 1, 1);
    stem(x, 'filled', 'LineStyle', 'none', 'MarkerSize', 8);
    title('Input x(t) of the Convolution');
    xlabel('Time (s)');
    ylabel('Magnitude');
    xline(0);             
    yline(0);             

    subplot(2, 1, 2);
    stem(y, 'filled', 'LineStyle', 'none', 'MarkerSize', 8);
    title('Output y(t) of the Convolution for h = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0]');
    xlabel('Time (s)');
    ylabel('Magnitude');
    xline(0);             
    yline(0);             
end

function define_plot_X_H3()
    x = [0, 1, 1, 1, 1, 0, 0, 0, 0, 0];
    h1 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    h2 = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
    h3 = [1, 0, 1, 0, 0, 0, 0, 0, 0, 0];

    y1 = myConv(x, h1);
    y2 = myConv(x, h2);
    y3 = myConv(x, h3);

    yT = y1 + y2;

    figure(Name='Exercise 2)D)');
    subplot(3, 1, 1);
    stem(x, 'filled', 'LineStyle', 'none', 'MarkerSize', 8);
    title('Input x(t) of the Convolution');
    xlabel('Time (s)');
    ylabel('Magnitude');
    xline(0);             
    yline(0);             

    subplot(3, 1, 2);
    stem(y3, 'filled', 'LineStyle', 'none', 'MarkerSize', 8);
    title('Output y(t) of the Convolution for h3=[1, 0, 1, 0, 0, 0, 0, 0, 0, 0]');
    xlabel('Time (s)');
    ylabel('Magnitude');
    xline(0);             
    yline(0);             

    % This point tries to make a point where the convolution for h3=[1,0,1,...] is equal to the convolution of h1=[1,0,0,...] + convolution of h2=[0,0,1,...]
    subplot(3, 1, 3);
    stem(yT, 'filled', 'LineStyle', 'none', 'MarkerSize', 8);
    title('Sum yT(t) of the convolution outputs h1=[1, 0, 0, 0...] + h2=[0, 0, 1, 0, ...]');
    xlabel('Time (s)');
    ylabel('Magnitude');
    xline(0);            
    yline(0);             
end