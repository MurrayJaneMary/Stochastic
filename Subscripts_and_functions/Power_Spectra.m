%Power spectra at CMB and surface using approach of Stacey (1992)

figure
fig = gcf; 
Re_Rcmb = (6371/3485);% Radius of Earth / Radius of Core

    MeanSqField_surface= zeros(max_degreePSV,size(gh_full,2));
    for degree = 1:max_degreePSV
        startrow = (degree+1)*(degree-1) + 1;
        endrow = (degree+2)*degree;
    
        for timestep = 1:Nt
            MeanSqField_surface(degree,timestep) = (degree+1)*sumsqr(gh_full(startrow:endrow,timestep));
            MeanSqField_CMB(degree,timestep) = MeanSqField_surface(degree,timestep) .* Re_Rcmb.^(2*degree+4);
    
        end
    end

    MeanSqField_TAF= mean(MeanSqField_surface,2);
    MeanSqField_TAF_CMB= mean(MeanSqField_CMB,2);
    
        semilogy(1:max_degreePSV, MeanSqField_TAF_CMB,'O-','LineWidth',2,...
            'MarkerEdgeColor','r','Color','r')
        hold on
        semilogy(1:max_degreePSV, MeanSqField_TAF,'O-','LineWidth',2,...
            'MarkerEdgeColor','b', 'Color','b')
        hold off
        xlabel('Degree'); ylabel('Power'); title('Power Spectra (Stacey, 1992)')
        MakePlotAxesNice