valid = Results.ADF_median ~= 0;
ADFratio = Results.ADF_iqr(valid) ./ Results.ADF_median(valid);
rev_rate = Results.RevRate(valid);

tiledlayout (2,3)

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
plot(rev_rate, Results.Model_G_a_S_nocut,  'go')
xlabel('RevRate (Myr^{-1})')
ylabel('Model_G_a')

title(' Model G A vs RevRate')

% Tile 5
nexttile
plot(rev_rate, Results.t_trans,  'mo')
xlabel('RevRate (Myr^{-1})')
ylabel('Transistion time fraction')

title(' Transistion time fraction')
