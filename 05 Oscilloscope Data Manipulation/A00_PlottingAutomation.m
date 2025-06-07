%% Purpose of the program: Import oscilloscope data with ease, using the raw file and skipping the first 20 lines that include configuration, and automating the plotting with desired settings.

%% One file example. Skip the first 20 rows so that row 21 becomes the header line, and read all data below it from a CSV file
% T = readtable("tek0000ALL.csv", ...
%               "HeaderLines", 20);

%% Multiple files
NUMBER_OF_FILES = 8;
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

% This program is configured to display information captured in XY mode.
% However, changing this to Y-t is extremely easy.
for i = 1:NUMBER_OF_FILES
    T = allData{i};
    x = T{:,2};
    y = T{:,3};
    
    txt = sprintf('(X-Y Mode) Oscilloscope Data from file %d', i);

    figure
    plot(x, y, '.-', 'LineWidth', 0.6, 'MarkerSize', 5);
    xlabel('Channel 1');
    ylabel('Channel 2');
    grid on;
    title(txt);
end