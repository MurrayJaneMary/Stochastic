valid = Results.ADF_median ~= 0;
ADFratio = Results.ADF_iqr(valid) ./ Results.ADF_median(valid);
rev_rate = Results.RevRate(valid);

tiledlayout (2,2)

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
plot(ADFratio, rev_rate, 'ko')
xlabel('ADF IQR / ADF Median')
ylabel('RevRate')
title('RevRate vs ADF IQR-to-Median Ratio')


% Tile 4
nexttile
plot(Results.Model_G_A, rev_rate, 'go')
xlabel('Model_G_A')
ylabel('RevRate')
title('RevRate vs Model G A')
