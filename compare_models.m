%Read in model metadata from Buffet model
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
    disp(modelName)
    filename = fullfile(modelName, [modelName, '_metadata.json']);

    if isfile(filename)
        jsonText = fileread(filename);
        metadata = jsondecode(jsonText);

        Amin(i) = metadata.Amin;
        Amax(i) = metadata.Amax;
        Recurrence(i) = metadata.Recurrence;
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


%Plot results 

valid = Results.ADF_median ~= 0;
ADFratio = Results.ADF_iqr(valid) ./ Results.ADF_median(valid);
rev_rate = Results.RevRate(valid);

valid = Results.VDM_median ~= 0;
VDMratio = Results.VDM_IQR(valid) ./ Results.VDM_median(valid);

tiledlayout (2,4)

% Tile 1
nexttile
plot(Results.ADF_median, Results.ADF_iqr, 'bo')
xlabel('ADF median')
ylabel('ADF IQR')
title('ADF IQR vs ADF Median')

% Tile 2
nexttile
plot(Results.RevRate, Results.ADF_median, 'ro')
xlabel('RevRate')
ylabel('ADF median')

title('ADF Median vs RevRate ')

% Tile 3
nexttile
plot(rev_rate, ADFratio,  'ko')
xlabel('RevRate (Myr^{-1})')
ylabel('ADF IQR / ADF Median')

title('ADF IQR-to-Median Ratio vs RevRate')


% Tile 4
nexttile
plot(rev_rate, Results.t_trans,  'mo')
xlabel('RevRate (Myr^{-1})')
ylabel('Transistion time fraction')

title(' Transistion time fraction')

% Tile 5
nexttile
plot(rev_rate, Results.Model_G_a_S_nocut,  'go')
xlabel('RevRate (Myr^{-1})')
ylabel('Model_G_a')

title(' Model G A vs RevRate')

% Tile 6
nexttile
plot(rev_rate, Results.VDM_median,  'mo')
xlabel('RevRate (Myr^{-1})')
ylabel('VDM median (Zam^{2})')
title(' VDM median')

% Tile 7
nexttile
plot(rev_rate, Results.VDM_IQR,  'mo')
xlabel('RevRate (Myr^{-1})')
ylabel('VDM IQR (Zam^{2})')
title(' VDM IQR')

%Tile 8
nexttile
plot(rev_rate, VDMratio, 'mo')
xlabel('RevRate (Myr^{-1})')
ylabel('VDM IQR/Median')
title(' VDM IQR to Median ratio vs Rev rate')