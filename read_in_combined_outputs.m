modelNames = {'model038a'};
for i = 1:length(modelNames)
    modelName = modelNames{i};
    filePath = fullfile('C:\Users\murray98\Documents\Bruce Buffett model\code', ...
                        modelName, 'combined_output_julia.csv');
    disp(filePath)
    metadata_filename = fullfile(modelName, [modelName, '_metadata.json']);

    if isfile(metadata_filename)
        jsonText = fileread(metadata_filename);
        metadata = jsondecode(jsonText);
    else
        warning('Missing metadata for %s', modelName);
    end
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
    timesteps = dataArray(1, :).';  % Transpose to column vector
    gh_full = dataArray(2:end, :);
    %gh_full = gh_full./1000 %conversion from nanoT to microT
    if isfield(metadata, 'scale') && metadata.scale == 1
        gh_full(1, :) = gh_full(1, :) * 5.70 / metadata.mean_x;
    end
    % Save timesteps and gh_full to separate .mat files in the same folder
    save(fullfile(folderPath, 'timesteps.mat'), 'timesteps');
    save(fullfile(folderPath, 'gh_full.mat'), 'gh_full');
end