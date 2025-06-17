% Extract your variables
x = FullResults(:,4)
y = FullResults(:,5)
% Define the model: A*sin(omega*x + phi) + C
ft = fittype('A*cos(w*x*pi/180 + phi) + C', ...
    'independent', 'x', ...
    'coefficients', {'A', 'w', 'phi', 'C'});


% Provide initial guesses for parameters to help convergence
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.StartPoint = [1, 1, 0, mean(y)];

% Fit the model
[fitresult, gof] = fit(x, y, ft, opts);

% Plot the fit
figure;
hold on;
plot(fitresult, x, y);

% Plot additional data
plot(PallatMOD2,SbprimeMOD2_Snc,'g-', 'DisplayName', 'x2 vs y2'); hold off


xlabel('Pallat');
ylabel('S (no cutoff)');
title('Sinusoidal Fit');
legend('Data', 'Sin curve', 'ModelG');
