n = height(Results);

Amin = zeros(n, 1);
Amax = zeros(n, 1);
Recurrence = zeros(n, 1);
MeanX = zeros(n, 1);
StdX = zeros(n, 1);
MedianX = zeros(n, 1);
IQRX = zeros(n, 1);
RevRate = zeros(n, 1);
PropReversing = zeros(n, 1);

for i = 1:n
    modelName = Results.Model{i};
    filename = fullfile(modelName, [modelName, '_metadata.json']);

    if isfile(filename)
        jsonText = fileread(filename);
        metadata = jsondecode(jsonText);

        Amin(i) = metadata.Amin;
        Amax(i) = metadata.Amax;
        Recurrence(i) = metadata.recurrence;
        MeanX(i) = metadata.mean_x;
        StdX(i) = metadata.std_x;
        MedianX(i) = metadata.median_x;
        IQRX(i) = metadata.iqr_x;
        RevRate(i) = metadata.rev_rate;
        PropReversing(i) = metadata.prop_reversing;
    else
        warning('Missing metadata for %s', modelName);
    end
end

Results.Amin = Amin;
Results.Amax = Amax;
Results.Recurrence = Recurrence;
Results.MeanX = MeanX;
Results.StdX = StdX;
Results.MedianX = MedianX;
Results.IQRX = IQRX;
Results.RevRate = RevRate;
Results.PropReversing = PropReversing;
