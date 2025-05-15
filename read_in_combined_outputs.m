% Define the path to the large CSV file
filePath = 'C:\Users\murray98\Documents\Bruce Buffett model\code\model3\combined_output_2025-05-13.csv';

% Get the folder of the input file to save output in the same location
[folderPath, ~, ~] = fileparts(filePath);

% Use a datastore to handle large number of columns
ds = tabularTextDatastore(filePath, ...
    'Delimiter', ',', ...
    'ReadVariableNames', false);  % Use false if you want to skip headers

% Read all the data
dataTable = readall(ds);

% Convert to numeric array (if all columns are numeric)
dataArray = table2array(dataTable);

% Extract the first row as timesteps and the rest as gh_full
timesteps = dataArray(1, :);
gh_full = dataArray(2:end, :);

% Save timesteps and gh_full to separate .mat files in the same folder
save(fullfile(folderPath, 'timesteps.mat'), 'timesteps');
save(fullfile(folderPath, 'gh_full.mat'), 'gh_full');
