% If the third parameter would be a cell instead of a table, it would be much easier to manipulate. Instead, we need to use the column names to operate it.
function newTable = uniformData(time_uni, time_s, myTable, widthTable)
    for i = 1 : widthTable
        name = myTable.Properties.VariableNames{i};
        originalColumn = myTable.(name);

        newColumn = interp1(...
        time_s,                  ...% original timestamps
        originalColumn(:),       ...% your array of data (DC removed)
        time_uni,                ...% uniform timestamps
        'linear');               ...% interpolation method

        newTable.(name) = newColumn;
    end
    
end