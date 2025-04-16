function [CentralLat, CentralLong, S, Sprime] = SphExpFunc2(VGPLat, VGPLong)

% Generates Spherical exponential central direction, Scatter (S) and S'(upscaled according to Doubrovine)
% uses fmincon to search for smallest sum of angular distance between
% central direction and vgps

%*****TESTED against Pavel's results in Test_of_SphExpFunc_WsCorr.xlsx

A = fisherfun([VGPLong, VGPLat]);
FishMeanLat = A(2);
FishMeanLong = A(1);
N = A(3);
% AngDist = zeros(N,1); %Preallocate memory
% AngDistSum = zeros(999999,3);
%Rotate VGPs such that fisher mean is at geog pole
AxisAz(1:N,1) = FishMeanLong+90;
AxisPl(1:N,1) = 0;
RotAng(1:N,1) = -90+FishMeanLat;

%[NewVGPLong, NewVGPLat] = SphereRot([VGPLong VGPLat AxisAz AxisPl RotAng]);

% new code by RKB
f = @(x)angdissum(x,VGPLat,VGPLong); % de-parameterize fmincon function, by adding individ. VGPs here, x is sum of ang. distance
options = optimset('TolFun', 1e-8, 'TolX', 1e-8, 'MaxIter', 5e3, 'Display', 'none');
% minimizes function f (sum of distance) with fisher mean as initial guess amd bounded over geographic range
[x,~] = fmincon(f,[FishMeanLat,FishMeanLong],[],[],[],[],[-90 0],[90 360],[],options); 

Best(1) = x(1); Best(2) = x(2);

% Best = AngDistSum2(size(AngDistSum2,1),:);

%[CentralLong, CentralLat] = SphereRot([Best(2) Best(1) AxisAz(1) AxisPl(1) -RotAng(1)]);

CentralLong = Best(2);
CentralLat = Best(1);

CentDir(1:N,1) = CentralLat; CentDir(1:N,2) = CentralLong;
AngDisSE = angdis([CentDir VGPLat VGPLong]);

%Calculate S using the SE central direction
AngDisSE_sqr = AngDisSE.^2;
SumAngDisSE_rms = (sum(AngDisSE_sqr));
S = sqrt((1/N).*SumAngDisSE_rms);

%Calculate S' using the SE central direction
SumAngDisSE = sum(AngDisSE);
Sprime_noUps = (1/N).*SumAngDisSE;
%Upscale by equation 12 from Doubrovine et al. (in prep)
Sprime = Sprime_noUps .* 1/(2-(N/(sqrt(N.*(N-1)))));

