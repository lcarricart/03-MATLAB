%% SS1 LAB2 Problem 1

%% User input 
prompt = {'Enter the value for t<0', 'Enter the value for t=0', 'Enter the value for t>0'};
dlgtitle = 'Signal Input';
dims = [1 35; 1 35; 1 35];
definput = {'0', '0.5', '1'};

input = inputdlg(prompt, dlgtitle, dims, definput);

%% Assigning the values
if isempty(input)
    t_negative_value = 0;
    t_0_value = 0.5;
    t_positive_value = 1;
else
    t_negative_value = str2double(input{1});
    t_0_value = str2double(input{2});
    t_positive_value = str2double(input{3});
end

%% The main function x(t)
N = 201;

t_negative = -(N-1)/2:0.01:0;
t_positive = 0:0.01:(N-1)/2;
t = [t_negative, NaN, 0, NaN, t_positive];

x_negative = ones(size(t_negative)) * t_negative_value;
x_positive = ones(size(t_positive)) * t_positive_value;
x = [x_negative, NaN, t_0_value, NaN, x_positive];

% Calculating the Odd and Even parts
% fliplr(x) reverses the vector from left to right.
xe = 0.5 * (x + fliplr(x));
xo = 0.5 * (x - fliplr(x));

%% Plotting
figure

subplot(2,2,1)
plot(t, x, '-', 0, t_0_value, 'o')
title('Original Signal: x(t)')
xlabel('Time (t)')
ylabel('x(t)')
ylim([-2, 2])
grid on

subplot(2,2,2)
plot(t, xe, '-', 0, xe(t == 0), 'o')
title('Even part: xe(t)')
xlabel('Time (t)')
ylabel('xe(t)')
ylim([-2, 2])
grid on

subplot(2,2,3)
plot(t, xo, '-', 0, xo(t == 0), 'o')
title('Odd part: xo(t)')
xlabel('Time (t)')
ylabel('xo(t)')
ylim([-2, 2])
grid on

subplot(2,2,4)
plot(t, xe + xo, '-', 0, t_0_value, 'o')
title('Sum: xe(t) + xo(t)')
xlabel('Time (t)')
ylabel('xe(t) + xo(t)')
ylim([-2, 2])
grid on