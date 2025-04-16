function [plat, plong] = polecalc(A)
    % POLECALC Calculate plate latitude and longitude from input data.
    %
    %   [plat, plong] = polecalc(A) calculates the plate latitude (plat) and
    %   longitude (plong) from input matrix A, where
    %     A(:,1) = declination in degrees,
    %     A(:,2) = inclination in degrees,
    %     A(:,3) = site lat in degrees,
    %     A(:,4) = site long in degrees.
    %
    %   Outputs:
    %     plat  - Plate latitude in degrees
    %     plong - Plate longitude in degrees

    decr   = deg2rad(A(:,1)); % Declination
    incr   = deg2rad(A(:,2)); % Inclination
    slatr  = deg2rad(A(:,3)); % Slatitude
    slongr = deg2rad(A(:,4)); % Slongitude

    % Calculate Magnetic Colatitude
    magcolat = pi/2 - atan(tan(incr) / 2);
    platr = asin(sin(slatr) .* cos(magcolat) + cos(slatr) .* sin(magcolat) .* cos(decr));

    % Ensure the argument of asin is within [-1, 1] to avoid complex numbers
    beta_arg = (sin(magcolat) .* sin(decr)) ./ cos(platr);
    beta_arg = max(min(beta_arg, 1), -1); % Clamp values to [-1, 1]
    beta = asin(beta_arg);

    maxrow = size(A, 1);
    plong = zeros(maxrow, 1);

    % near/farsided pole
    condition = cos(magcolat) >= (sin(slatr) .* sin(platr));

    plong(condition) = rad2deg(slongr(condition) + beta(condition));
    plong(~condition) = rad2deg(slongr(~condition) - beta(~condition)) + 180;

    % Ensure Longitude is Within [0, 360)
    plong = mod(plong, 360);

    plat = rad2deg(platr);
end
