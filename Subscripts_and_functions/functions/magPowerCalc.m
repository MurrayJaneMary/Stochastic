%Returns magnetic field power (mean square field) at depth given due to all terms per degree l using
%equations 24.21 and 24.23 of Stacey 1992
% send CMB=1 for field at CMB; otherwise will calc at surface

function [MeanSqField, MeanField] = magPowerCalc(gh_full, CMB, maxdegree)

Re_Rcmb = (6371/3485);% Radius of Earth / Radius of Core
MeanSqField = zeros(maxdegree,size(gh_full,2)); MeanSqField=[];


for degree = 1:maxdegree
    startrow = (degree+1)*(degree-1) + 1;
    endrow = (degree+2)*degree;

    for timestep = 1:size(gh_full,2)
        %For each timestep
    %Calculate mean square surface field due to all terms of degree l using
    %equation 24.21 of Stacey 1992
        
        MeanSqField(degree,timestep) = (degree+1)*sumsqr(gh_full(startrow:endrow,timestep));

        if CMB==1
        %Downward continue using equation 24.23
            MeanSqField(degree,timestep) = MeanSqField(degree,timestep) .* Re_Rcmb.^(2*degree+4);
        end
            MeanField(degree,timestep) = sqrt(MeanSqField(degree,timestep));
        
    end
end
