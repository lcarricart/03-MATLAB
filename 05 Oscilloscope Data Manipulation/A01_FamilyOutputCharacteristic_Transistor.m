%% Purpose of the program: Import data from the oscilloscope in X-Y mode and graph the family of output characteristics of a transistor.
% This data manipulation corresponds to Electronics 1, Laboratory 3, Exercise 1.

%% Multiple files
NUMBER_OF_FILES = 8;
colors = lines(NUMBER_OF_FILES);

% Clears out the screen and figures
clc;
close all;

% Preallocate cell array to hold all imported files
allData = cell(8,1);

% Iterate the reading process
for i = 0:(NUMBER_OF_FILES - 1)
    % Build filename with zero-padded index
    fname = sprintf("tek%04dALL.csv", i);
    
    % Read table, skipping the first 20 rows so row 21 is header
    T = readtable(fname, "HeaderLines", 20);
    
    % Store in cell array
    allData{i+1} = T;
    
    % % (Optional) Display first 5×5 block to verify
    % fprintf("First 5×5 of %s:\n", fname);
    % disp(T(1:min(5,height(T)), 1:min(5,width(T))));
end

figure;
for i = 1:NUMBER_OF_FILES
    T = allData{i};
    x = T{:,2};
    y = T{:,3};

    txt = sprintf('Family of Output Characteristics');

    plot(x, y, '.-', 'Color', colors(i,:), 'LineWidth',0.6, 'MarkerSize',5);
    hold on;
end

xlabel('U_{CE} []');
ylabel('I{C} []');
hold off;
title(txt);