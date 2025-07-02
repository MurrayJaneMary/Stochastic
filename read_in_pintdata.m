
dataDir = 'C:\Users\murray98\Documents\Bruce Buffett model\code';
files = dir(fullfile(dataDir, 'PINTdata*.json'));

for k = 1:length(files)
    filename = fullfile(dataDir, files(k).name);
    disp(['Reading in: ', filename])

    jsonText = fileread(filename);
    metadata = jsondecode(jsonText);

    % Convert to one-row table
    T = struct2table(metadata, 'AsArray', true);

    % Align variable names: add missing columns to both sides
    allVars = union(Results.Properties.VariableNames, ...
                    T.Properties.VariableNames);

    % For missing columns in Results, fill with default values for all rows
    for v = setdiff(allVars, Results.Properties.VariableNames)
        Results.(v{1}) = repmat(missingValue(T.(v{1})), height(Results), 1);
    end

    % For missing columns in T, fill with default values for all rows (usually 1 row)
    for v = setdiff(allVars, T.Properties.VariableNames)
        T.(v{1}) = repmat(missingValue(Results.(v{1})), height(T), 1);
    end

    % Reorder T to match Results
    T = T(:, Results.Properties.VariableNames);

    % Append
    Results = [Results; T];
end

function val = missingValue(col)
    if iscell(col)
        val = {[]};
    elseif isstring(col) || ischar(col)
        val = "";
    elseif isnumeric(col)
        val = NaN;
    elseif islogical(col)
        val = false;
    else
        val = missing;
    end
end

