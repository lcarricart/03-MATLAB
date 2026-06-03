clc; clear;
 function y = myConv (x, h)
 % output length
 Nx = length(x);
 Nh = length(h);
 Ny = Nx + Nh - 1;
 
 % Before looping, creating a zero vector
 y = zeros(1, Ny);

 for n=1:Ny
     for k =1:Nx
         h_index = n-k+1;
         if h_index >=1 && h_index<=Nh
             y(n) = y(n) + x(k) * h(h_index);
         end
     end
 end
end

% Defining the input signal x (rectangular pulse)
x = zeros(1, 10);
x(2:5) = 1;    % indices 2,3,4,5 are 1

% Defining the impulse response h
h = zeros(1, 10);
h(1) = 1;        % unit impulse at index 1

% Perform convolution using function from Part (a)
y = myConv(x, h);

%part c
h = zeros(1,10);
h(3) = 1;

y2 = myConv(x,h);

% part d
h = zeros(1,10);
h(1) = 1;
h(3) = 1;

y3 = myConv(x, h);

% Ploting the input and output
figure;

subplot(4,1,1);
stem(x, 'filled');
title('Input Signal x[n]');
xlabel('n');
ylabel('Amplitude');

subplot(4,1,2);
stem(y, 'filled');
title('Output Signal y[n] = x[n] * h[n]');
xlabel('n');
ylabel('Amplitude');

subplot(4,1,3);
stem(y2, 'filled');
title('Part (c): h(3)=1 -> y is shifted copy of x');
xlabel('n');
ylabel('Amplitude');


subplot(4,1,4);
stem(y3, 'filled');
title('Part (d): y3 = x[n] + x[n-2] (two impulses)');
xlabel('n');
ylabel('Amplitude');

