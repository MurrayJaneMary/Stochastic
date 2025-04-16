function PSVindex = PSVindexCalc(VGPlat, VDM, M0)

%Calculates PSVindex from Panovska et al. (2017) using default value of
%M0=80 ZAm^2 if nothing is passed.

if exist("M0","var") == 0
    M0 = 80;
end

if VGPlat > 0
    Theta_prime = 90-VGPlat;
else
    Theta_prime = 90+VGPlat;
end

PSVindex = (Theta_prime .* M0)/(180 .* VDM);
