function DipoleMoment = DipoleMomentFunc(g10,g11,h11)

%function to return dipole moment (in ZAm^2) from relevant gauss coefficients in nanoTesla.
% for ADM, leave last two input terms as zero

g10=g10 ./ 10^9; 
g11=g11 ./ 10^9; 
h11=h11 ./ 10^9; %convert to Tesla

Radius = 6370e3;
mu_0 = 1.257e-6;
first = (4*pi*Radius^3)/mu_0;
second = sqrt(g10.^2 + g11.^2 + h11.^2);

DipoleMoment = first .* second ./ 10^21;