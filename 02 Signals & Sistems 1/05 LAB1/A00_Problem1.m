%% Purpose of the program: SS1, Laboratory 1.

%% Problem 1: "Create the following vectors"
vectorA = linspace(1, 5, 5); % from 1 to 5, 5 values
vectorB = [3, -1, 2]';       % transpose a vector
vectorC = 0:0.1:10;          % from 0 to 10 in steps of 0.1
vectorD = 10:-0.1:0;         % opposite, mind the -1!
disp('Problem 1: vectors created')

% Test! Display the created vectors
% disp('Vector A:'); disp(vectorA);
% disp('Vector B:'); disp(vectorB);
% disp('Vector C:'); disp(vectorC);
% disp('Vector D:'); disp(vectorD);

%% Problem 2: "Use the Matlab function acos to calculate the angle between the vectors (4,2,5) and (3,8,9).
vector1 = [4, 2, 5];
vector2 = [3, 8, 9];
angleVector = acos(dot(vector1, vector2) / ...
                   (norm(vector1) * norm(vector2)));

% Display the calculated angle in radians
disp('Problem 2');
disp('Angle between the vectors (in radians):'); disp(angleVector);
disp('Angle between the vectors (in degrees):'); disp(rad2deg(angleVector));

%% Problem 3: Calculus of sums
%syms n k
%result = symsum(n^k, k, 1, 10);
disp('Problem 3');
%disp('[Example] The result of the sum "n^k" from 1 < k < 10 is:'); disp(result);

%sumA = symsum(k, k, 0, 100);
%sumB = symsum(k^2, k, 0, 100);
%disp('Sum A (k) = '); disp(sumA);
%disp('Sum B (k) = '); disp(sumB);

%% Problem 4:
k = (1:20).';
sineTaylor = @(x) sum( (-1).^(k-1) .* x.^(2*k-1) ./ factorial(2*k-1), 1);
x = linspace(-4*pi, 4*pi, 3000);
figure(1);
plot(x, sineTaylor(x));
xlim([-5, 15]);
title('Taylor expansion of the sine function');
grid on;
disp('Problem 4');
disp('sineTaylor(-pi/2) = '); disp(sineTaylor(-pi/2));

% Adjusting to the exercise's requirements. "Evaluate this series for x = -pi/2 by summing up the terms until the magnitude of the last summand drops below the value 10^(-8)"
value = -pi/2;
functionProblem4 = 0;
lastSummand = 1;                        % Arbitrary initialization (since I would need a do-while)
i = 1;

while (abs(lastSummand) >= 10^(-8))     % The abs() function is very important because lastSummand goes + and -
    lastSummand = ((-1)^(i-1) * value^(2*i-1)) / factorial(2*i-1);
    functionProblem4 = functionProblem4 + lastSummand;
    i = i + 1;
end

disp('The result of Problem 4 by the means stated in the exercise is:');
disp(functionProblem4);

%% Problem 5: "Write a function mysum that takes arguments for two numbers and returns their sum"
function [result] = mySum(a, b)
    result = a + b;
end

resultProblem5 = mySum(5, 4);
disp('Problem 5:'); disp('mySum(5, 4) --> 5 + 4');
disp(resultProblem5);

%% Problem 6: "Plot the following two signals in the time range between t = 0s and t = 2s"
t = 0:0.01:2;
freq1 = 2; % Hz
freq2 = 5; % Hz

function1_P6 = @(t) cos(2*pi*freq1*t);
%function2_P6 = @(t) (t<1).*cos(2*pi*freq1.*t) + (t>=1).*cos(2*pi*freq2.*t);
function2_P6 = @(t) (t<1).*cos(2*pi*freq1.*t);
function3_P6 = @(t) (t>=1).*cos(2*pi*freq2.*t);
function4_P6 = [function2_P6, function3_P6];

figure(2);
plot(t, function1_P6(t), 'r', t, function4_P6(t), 'b');
title('Problem 6 Plotting');
grid on;

