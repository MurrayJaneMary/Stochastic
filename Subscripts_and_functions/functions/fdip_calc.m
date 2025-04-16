function [fdip,fdip_timesteps] = fdip_calc(gh_full,maxdegree)
%Correct - tested

%load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Geodynamo_Sims\nm.mat')

%Definition from Christensen and Aubert
% we determine the relative
% dipole field strength f dip, defined as the time-average ratio on the
% outer shell boundary of the mean dipole field strength to the field
% strength in harmonic degrees  = 1–12,


%Define maxdegree
%maxdegree = 10;

Re_Rcmb = (6371/3485);% Radius of Earth / Radius of Core
maxrow = (maxdegree+2)*maxdegree;

%For each timestep
%Calculate mean square surface field due to all terms of degree l using
%equation 24.21 of Stacey 1992




MeanSqField_surface = zeros(maxdegree,size(gh_full,2)); MeanSqField_surface=[];
MeanSqField_CMB = zeros(maxdegree,size(gh_full,2)); MeanSqField_CMB=[];
Total_Field = zeros(1,size(gh_full,2)); Total_Field = [];


for degree = 1:maxdegree
    startrow = (degree+1)*(degree-1) + 1;
    endrow = (degree+2)*degree;

    for timestep = 1:size(gh_full,2)
        %For each timestep
    %Calculate mean square surface field due to all terms of degree l using
    %equation 24.21 of Stacey 1992
        
        MeanSqField_surface(degree,timestep) = (degree+1)*sumsqr(gh_full(startrow:endrow,timestep));

        
        %Downward continue using equation 24.23
        MeanSqField_CMB(degree,timestep) = MeanSqField_surface(degree,timestep) .* Re_Rcmb.^(2*degree+4);
        
        %Take sqrt to get mean field strength
        MeanField_CMB(degree,timestep) = sqrt(MeanSqField_CMB(degree,timestep));

        %Sum from degree 1 to maxdegree - this is total magnetic field strength
        Total_Field(1,timestep) = sqrt(sumsqr(MeanField_CMB(:,timestep)));

    end
end


%Divide mean field strengh of l=1 terms by total field strength to obtain
%relative dipole strength
fdip_timesteps = MeanField_CMB(1,:) ./ Total_Field(1,:);

%Average across all timesteps to get f_dip
fdip = mean(fdip_timesteps);

