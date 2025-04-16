function [timestep, slat, slong, dec, inc, F, plat, plong, VDM, VADM]=QueryModel(gh_full,timestep,slat,slong,max_degree)

%generic query of gh at a time and place. 
%Use convert_gh.m to generate gh_full

gh = gh_full(:,timestep); %column vector of Gauss coeffs at that timestep
nlat=slat; elong=slong; alt=0;

 

B = g2fxFunc_vRKB(gh,max_degree,nlat,elong,alt); %calls script to get B

Bcart = Cart2DI(B');
    dec = Bcart(1,1); inc = Bcart(1,2); F=Bcart(1,3).*1000000; %gives F in micro_T
    
    polecalcin =[dec, inc, slat, slong];
    [plat, plong] = polecalc(polecalcin);
    
    VDM = VDMcalc(F,inc); %in ZAm^2
    inclat=lat2inc(nlat);
    VADM = VDMcalc(F,inclat); %in ZAm^2
   
    %output(row,:) = [i_VDM nlat elong dec inc F VDM VADM];