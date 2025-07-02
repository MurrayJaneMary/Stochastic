%Plot results 
% Separate indices by modelname prefix
isPINT = startsWith(string(Results.Model), "PINTdata");
isModel = startsWith(string(Results.Model), "model");

valid = Results.ADF_median ~= 0;
ADFratio = Results.ADF_iqr(valid) ./ Results.ADF_median(valid);
rev_rate = Results.RevRate(valid);

valid = Results.VDM_median ~= 0;
VDMratio_PINT = Results.VDM_IQR(isPINT) ./ Results.VDM_median(isPINT);
VDMratio_model = Results.VDM_IQR(isModel) ./ Results.VDM_median(isModel);


tiledlayout (2,4)



% Tile 1
nexttile
plot(Results.ADF_median, Results.ADF_iqr, 'bo')
xlabel('ADF median')
ylabel('ADF IQR')
title('ADF IQR vs ADF Median')

% Tile 2
nexttile
plot(Results.RevRate, Results.ADF_median, 'bo')
xlabel('RevRate')
ylabel('ADF median')

title('ADF Median vs RevRate ')

% Tile 3
nexttile
plot(rev_rate, ADFratio,  'bo')
plot(Results.VDM_median(isPINT), Results.Model_G_a_S_nocut(isPINT), 'rs')
plot(Results.VDM_median(isModel), Results.Model_G_a_S_nocut(isModel), 'bo')
hold off
xlabel('VDM median (Zam^{2})')
ylabel('Model G A')
title('Model G vs VDM Median')


% Tile 4
nexttile
plot(rev_rate, Results.t_trans,  'bo')
xlabel('RevRate (Myr^{-1})')
ylabel('Transistion time fraction')

title(' Transistion time fraction')

% Tile 5
nexttile
hold on
plot(Results.RevRate(isPINT), Results.Model_G_a_S_nocut(isPINT), 'rs')
plot(Results.RevRate(isModel), Results.Model_G_a_S_nocut(isModel), 'bo')
hold off
xlabel('RevRate (Myr^{-1})')
ylabel('Model_G_a')

title(' Model G A vs RevRate')

% Tile 6
nexttile
hold on
plot(Results.RevRate(isPINT), Results.VDM_median(isPINT), 'rs')
plot(Results.RevRate(isModel), Results.VDM_median(isModel), 'bo')
hold off

xlabel('RevRate (Myr^{-1})')
ylabel('VDM median (Zam^{2})')
title(' VDM median')

% Tile 7
nexttile
hold on
plot(Results.RevRate(isPINT), Results.VDM_IQR(isPINT), 'rs')
plot(Results.RevRate(isModel), Results.VDM_IQR(isModel), 'bo')
hold off
xlabel('RevRate (Myr^{-1})')
ylabel('VDM IQR (Zam^{2})')
title(' VDM IQR')

%Tile 8
nexttile
hold on
plot(Results.RevRate(isPINT), VDMratio_PINT, 'rs', 'DisplayName', 'PINTdata')
plot(Results.RevRate(isModel), VDMratio_model, 'bo', 'DisplayName', 'model')
hold off
xlabel('RevRate (Myr^{-1})')
ylabel('VDM IQR / VDM Median')
title('VDM IQR to Median Ratio vs RevRate')
legend('Location', 'best')