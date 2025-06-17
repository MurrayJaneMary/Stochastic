function [ADF,NADF,Total_Field] = NADF_calc(gh_full,maxdegree)
%Calculate rms Axial Dipole field and rms non-axial dipole field at Earth's surface
%input GC's in nT and obtain fields in nT

%DF and NDF tested using IGRF2005 and values in Stacey 1992. ADF and NADF
%Ratio of ADF.^2/NADF.^2 = AD/NAD in Biggin et al. (2020)

%load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Geodynamo_Sims\nm.mat')

%Definition from Christensen and Aubert
% we determine the relative
% dipole field strength f dip, defined as the time-average ratio on the
% outer shell boundary of the mean dipole field strength to the field
% strength in harmonic degrees  = 1–12,


%Define maxdegree
maxdegree = 10;

Re_Rcmb = (6371/3485);% Radius of Earth / Radius of Core
maxrow = (maxdegree+2)*maxdegree;

%For each timestep
%Calculate mean square surface field due to all terms of degree l using
%equation 24.21 of Stacey 1992


MeanSqField_surface = zeros(maxdegree,size(gh_full,2)); MeanSqField_surface=[];
MeanField_surface = zeros(maxdegree,size(gh_full,2)); MeanField_surface=[];
%MeanSqField_CMB = zeros(maxdegree,size(gh_full,2)); MeanSqField_CMB=[];

% Total_Field = zeros(1,size(gh_full,2)); Total_Field = [];
% Total_Field_NADF = zeros(1,size(gh_full,2)); Total_Field_NADF = [];

gh_full_NADF = gh_full;
gh_full_NADF(1,:) = zeros(1,size(gh_full,2));


for degree = 1:maxdegree
    startrow = (degree+1)*(degree-1) + 1;
    endrow = (degree+2)*degree;

    for timestep = 1:size(gh_full,2)
                
        MeanSqField_surface = (degree+1)*sumsqr(gh_full(startrow:endrow,timestep));
        MeanSqField_surface_NADF = (degree+1)*sumsqr(gh_full_NADF(startrow:endrow,timestep));
        %Take sqrt to get mean field strength
        MeanField_surface(degree,timestep) = sqrt(MeanSqField_surface);
        MeanField_surface_NADF(degree,timestep) = sqrt(MeanSqField_surface_NADF);

        %Downward continue using equation 24.23
        %MeanSqField_CMB(degree,timestep) = MeanSqField_surface(degree,timestep) .* Re_Rcmb.^(2*degree+4);
        %MeanField_CMB(degree,timestep) = sqrt(MeanSqField_CMB(degree,timestep));
        
    end
end

Total_Field = sqrt(sum(MeanField_surface.^2));
DF = MeanField_surface(1,:);
NDF = sqrt(sum(MeanField_surface(2:10,:).^2));
%NADF = sqrt(sum(MeanField_surface_NADF.^2));
ADF = abs(sqrt(2).* gh_full(1,:));

NADF = sqrt(Total_Field.^2 - ADF .^2); %gives same answer as above but easier to write


% %Calculate total NADF  at each timestep
% for timestep = 1:size(gh_full,2)
%     %Sum from degree 1 to maxdegree - this is total magnetic field strength
%     Total_Field(1,timestep) = sqrt(sumsqr(MeanField_CMB(:,timestep)));
% 
%     g10_MeanSqField_surface(timestep) = (1+1)*sumsqr(gh_full(1,timestep));
%     g10_MeanSqField_CMB(timestep) = g10_MeanSqField_surface(timestep) .* Re_Rcmb.^(2*1+4);
%     g10_MeanField_CMB(timestep) = sqrt(g10_MeanSqField_CMB(timestep));
% end
% 
%     %Divide mean field strengh of l=1 terms by total field strength to obtain
% %relative dipole strength
% fdip_timesteps = MeanField_CMB(1,:) ./ Total_Field(1,:);
% %repeat for just axial dipole
% faxdip_timesteps = g10_MeanField_CMB(1,:) ./ Total_Field(1,:);
% 
% %Average across all timesteps to get f_dip
% fdip = mean(fdip_timesteps);
% 
% 
% %Average across all timesteps to get faxdip
% faxdip = mean(faxdip_timesteps);
