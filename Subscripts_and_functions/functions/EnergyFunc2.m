function [AxDipSurf, AxDip, AD_NADSurf, AD_NAD, Odd_EvenSurf, Odd_Even]= EnergyFunc2(gh,n,m);%, Odd_Evensurf, AD_NAD, AD_NADsurf, Odd_ADEvensurf, Even_ADOddsurf, AD_Oddsurf, AD_Evensurf, AD_Allsurf, Odd_Allsurf, Even_Allsurf] = EnergyFunc(gh,n,m)
%Function to calculate the energy of terms at the surface and CMB and then
%various ratios
% Following approach of Coe and Glatzmaier (2006) but to n,m =10

%This function is just for single timestep. For full time series with
%averages etc see EnergyFunc.m




Re_Rcmb = 6371/3485; %Radius of Earth/ Radius of Core

ii = 1; jj=1; NMgh = nan(20,3); W=[]; WSurf=[];


%find energy at CMB and surface
while ii < 48%120
    
   if m(ii) == 0 %axial component so no h component
       %Wsurf(jj,:) = gh(ii,:).^2; 
       WSurf(jj,:) = (n(ii)+1) .* gh(ii,:).^2;
       W(jj,:) = WSurf(jj,:) .* (Re_Rcmb.^(2.*(n(ii)+2)));
       NMgh(jj,1) = n(ii); NMgh(jj,2) = m(ii); %Record degree and order
       ii = ii+1;  jj=jj+1;
   else
       %Wsurf(jj,:) = gh(ii,:).^2 + gh(ii+1,:).^2;
       WSurf(jj,:) = (n(ii)+1) .* (gh(ii,:).^2 + gh(ii+1,:).^2);
       W(jj,:) = WSurf(jj,:) .* (Re_Rcmb.^(2.*(n(ii)+2)));
       NMgh(jj,1) = n(ii); NMgh(jj,2) = m(ii); %Record degree and order
       ii = ii+2; jj=jj+1;
   end
   
      
end

%Calculate time-instaneous ratios
AxDipSurf_num = sqrt(WSurf(1,:));
AxDip_num =sqrt(W(1,:));
AxDipSurf_denom = sqrt(sum(WSurf(:,:)));
AxDip_denom = sqrt(sum(W(:,:)));

AxDipSurf = AxDipSurf_num ./ AxDipSurf_denom;
AxDip = AxDip_num ./ AxDip_denom;

WSurf_total = sum(WSurf,1);
W_total = sum(W,1);

AD_NADSurf(1,:) = WSurf(1,:) ./ (WSurf_total(1,:) - WSurf(1,:));
AD_NAD(1,:) = W(1,:) ./ (W_total(1,:) - W(1,:));

%Calculate median ratios
AxDipSurf_Med = median(AxDipSurf);
AxDip_Med = median(AxDip);
AD_NADSurf_Med = median(AD_NADSurf);
AD_NAD_Med = median(AD_NAD); 



%EvenW = zeros(65,1); OddW = zeros(65,1); EvenWsurf = zeros(65,1); OddWsurf = zeros(65,1);
EvenW = 0; OddW = 0; EvenWsurf = 0; OddWsurf = 0;
for kk = 2:27%65
%for kk = 1:65
    SumMN = NMgh(kk,1)+NMgh(kk,2);
    if floor(SumMN./2) == SumMN./2 %check even
        
        EvenW(1,1) = EvenW(1,1) + W(kk,1);
        EvenWsurf(1,1) = EvenWsurf(1,1) + WSurf(kk,1);
    else
        OddW(1,1) = OddW(1,1) + W(kk,1);
        OddWsurf(1,1) = OddWsurf(1,1) + WSurf(kk,1);
    end
end

Odd_Even = OddW./EvenW;
Odd_EvenSurf = OddWsurf./ EvenWsurf;

Odd_EvenSurf_Med  = median(Odd_EvenSurf);
Odd_Even_Med = median(Odd_Even);

