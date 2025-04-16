% Define the path to the large CSV file
filePath = 'C:\Users\murray98\Documents\Bruce Buffett model\code\model2\combined_output.csv';

% Use a datastore to handle large number of columns
ds = tabularTextDatastore(filePath, ...
    'Delimiter', ',', ...
    'ReadVariableNames', false);  % Use false if you want to skip headers

% Preview or check variable names
% preview(ds)

% Read all the data
dataTable = readall(ds);

% Convert to numeric array (if all columns are numeric)
dataArray = table2array(dataTable);

% Save the array to a .mat file
save('combined_output.mat', 'dataArray');