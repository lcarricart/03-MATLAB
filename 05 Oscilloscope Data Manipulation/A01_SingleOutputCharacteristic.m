%% Purpose of the program: Import data from the oscilloscope in X-Y mode and graph ONE single set of the family of output characteristics of a transistor.
% This data manipulation corresponds to Electronics 1, Laboratory 3, Exercise 1.

%% Multiple files
NUMBER_OF_FILES = 1;
colors = lines(NUMBER_OF_FILES);
R  = 10e3;                 % 10 kΩ

% Clears out the screen and figures
clc;
%close all;

% Preallocate cell array to hold all imported files
allData = cell(1,1);

% Build filename with zero-padded index
fname = sprintf("tek%04dALL.csv", 1);

% Read table, skipping the first 20 rows so row 21 is header
T = readtable(fname, "HeaderLines", 20);

% Store in cell array
allData{1} = T;

% % (Optional) Display first 5×5 block to verify
% fprintf("First 5×5 of %s:\n", fname);
% disp(T(1:min(5,height(T)), 1:min(5,width(T))));

figure;
T = allData{1};
x = T{:,2};
y = T{:,3};

% In this laboratory, Y was measured to be the voltage over a resistor, to later calculate its current. I = V / R
y  = y / R;              
y_uA = y * 1e6;            

txt = sprintf('');

plot(x, y_uA, '.-', 'Color', colors(1,:), 'LineWidth',0.6, 'MarkerSize',3);
hold on;

xlabel('U_{CE} [V]');
ylabel('I{C} [uA]');
hold off;
title(txt);