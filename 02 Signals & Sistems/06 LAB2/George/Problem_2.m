%% SS1 LAB2 Problem 2

%% User input 
prompt = {'Enter disrcete time value', 'Enter frequency of the function in Hz', 'Enter the number of coefficients to be computed'};
dlgtitle = 'Fourier Series Input';
dims = [1 35; 1 35; 1 35];
definput = {'5', '1', '1000'};

input = inputdlg(prompt, dlgtitle, dims, definput);


%% Assigning the values
if isempty(input)
    time_lim = 5;
    f0 = 1;
    N = 1000;
else
    time_lim = str2double(input{1});
    f0 = str2double(input{2});
    N = str2double(input{3});
end

%% Values initialization
% a_k = 0;

num_points = max(5001, 1000 * time_lim * f0); % ensures that at least 5001 points will be plotted
t = linspace(0, time_lim, num_points);        % creating a vector
x_t = zeros(size(t));
w0 = 2*pi*f0;

%% Calculation of Fourier Series
for k = 1:2:N
    b_k = 4 / (pi * k);
    x_t = x_t + b_k * sin(w0 * k * t);
end

%% Plotting
plot(t, x_t)
xlabel('Time (s)')
ylabel('X(t)')
ylim([-1.5 1.5]);
title('Fourier Series representation of Rectangular Signal')
