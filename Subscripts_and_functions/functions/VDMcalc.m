%Calculates VDMs in ZAm^2 from F in micro_T and inc in degrees

function V = VDMcalc(F,I)

%perm of free space
mu = 4*pi()*1e-7;
a = 4*pi()*6371e3^3 / mu;

colat = 90 - degrees(atan(tan(radians(I)/2)));

b = sqrt(1+3.*(cos(radians(colat)).^2));

ab = a./b;

V = (ab .* F .*1e-6 .*1e-21); % assuming F is in microT, V is in ZAm2