% Plot:
% 1. VADM timeseries and median ADM
% 2. VADM distribution and iqr
% 3. VGP latitude timeseries; median VGP colat (polarity normalised);
% 4. VGP latitude distibution; iqr VGP colat (polarity normalised);
% 5. PSVindex time series; median PSVindex
% 6. PSVindex time series; iqr PSVindex


str1=sprintf('%0.0f',Ref_SiteLat);
str2=sprintf('%0.0f',Ref_SiteLong);
PlotFileName = strcat("Global_TimeSeries_Plots_Lat",str1,"_Long",str2,".fig");

figure;
%set(gcf, 'Position',  [100, 100, 500, 250.*Nmodels]) %set final fig to 250 x Nmodels
fig = gcf;
fig.WindowState = 'maximized';
t=tiledlayout(3,2);

%%

Nt = size(timesteps_dim,1);

gh_full_ADM = gh_full;% .* gh_CorrFac_DM;

for ts = 1: Nt
    
    gh = gh_full_ADM(:,ts);
    
    [dec(ts), inc(ts), F(ts), plat(ts), plong(ts), VDM(ts), VADM(ts)]=...
        QueryModel2(gh_full_ADM,ts,Ref_SiteLat,Ref_SiteLong,max_degreePSV);
    
        
end

VDM=VDM .* gh_CorrFac_VDM;
VADM=VADM .* gh_CorrFac_VDM;
%%

%1. VADM timeseries and median ADM
VADM_median = median(VADM);
VADM_median_line = [VADM_median VADM_median];

ax1= nexttile; %xlim([0,max(timesteps_dim)]);
plot(ax1,timesteps_dim,VADM);hold on;
plot(ax1, [timesteps_dim(1,1), timesteps_dim(Nt,1)] , VADM_median_line);hold off;
xlim([min(timesteps_dim),max(timesteps_dim)]);

xlabel(Time_axis_label);
if timesteps_dim(1) > timesteps_dim(2)
    plot(ax1,flip(timesteps_dim),flip(VADM));hold on;
    plot(ax1, flip([timesteps_dim(1,1), timesteps_dim(Nt,1)]) , VADM_median_line);hold off;
    set(gca, 'XDir', 'reverse'); xlabel(Time_axis_label);
end

ylim([0 1.1*max(VADM)]);
ylabel('Virtual Axial Dipole Moment (ZAm^2)');


str=sprintf('%2.0f',VADM_median);
Header = strcat('VADM_m_e_d_i_a_n = ',str);
title(Header);
MakePlotAxesNice

%%
% 2. VADM distribution and iqr
VADM_iqr = iqr(VADM);

ax1 = nexttile;
hist(ax1,VADM,100)
xlim([0,1.1*max(VADM)]);

str=sprintf('%2.0f',VADM_iqr);
Header = strcat('VADM_i_q_r = ',str);
title(Header);
xlabel('Virtual Axial Dipole Moment (ZAm^2)'); ylabel('N')
MakePlotAxesNice
%%
% 3. VGP lat timeseries

[p_colat] = 90-abs(plat); 
VGP_colat_median = median(p_colat);
VGP_lat_median_line = [90-VGP_colat_median 90-VGP_colat_median];

ax1= nexttile;
plot(ax1, timesteps_dim,plat,'Color','k')
xlim([min(timesteps_dim),max(timesteps_dim)]); ylim([-90 90]);yticks(-90:45:90);
hold on
plot(ax1, [timesteps_dim(1,1), timesteps_dim(Nt,1)] , VGP_lat_median_line);
hold off
xlabel(Time_axis_label);

if timesteps_dim(1) > timesteps_dim(2)

    plot(ax1,flip(timesteps_dim),flip(plat));hold on;
    plot(ax1, flip([timesteps_dim(1,1), timesteps_dim(Nt,1)]) , VGP_lat_median_line);hold off;
    set(gca, 'XDir', 'reverse'); xlabel(Time_axis_label);
end

ylabel('VGP latitude (°)');
%title('Dipole Tilt')
str=sprintf('%0.2f',VGP_colat_median);
Header = strcat('VGP Colat_m_e_d_i_a_n = ',str);
title(Header);
MakePlotAxesNice
%%
% 4. VGP lat distibution; iqr VGP colat
VGP_Colat_iqr = iqr(p_colat);

ax1 = nexttile;
hist(ax1,p_colat,100)

str=sprintf('%0.2f',VGP_Colat_iqr);
Header = strcat('VGP Colat_i_q_r = ',str);
title(Header);xlabel('VGP Colatitude (°)'); ylabel('N')
MakePlotAxesNice

%%



% 5. PSVindex time series; median PSVindex

for ts = 1:Nt
    PSVindex(ts) = PSVindexCalc(plat(ts), VADM(ts), VADM_median);
end

PSVi_median = median(PSVindex);
PSVi_median_line = [PSVi_median PSVi_median];


ax1= nexttile;
plot(ax1, timesteps_dim,PSVindex,'Color','k')
xlim([min(timesteps_dim),max(timesteps_dim)]); 
ylim([0 1.1*max(PSVindex)]);
ylabel('PSVindex_M_a_v_g');


hold on
plot(ax1, [timesteps_dim(1,1), timesteps_dim(Nt,1)] , PSVi_median_line);
hold off
xlabel('Time (Myr)');

if timesteps_dim(1) > timesteps_dim(2)

    plot(ax1,flip(timesteps_dim),flip(plat));hold on;
    plot(ax1, flip([timesteps_dim(1,1), timesteps_dim(Nt,1)]) , PSVi_median_line);hold off;
    set(gca, 'XDir', 'reverse'); xlabel(Time_axis_label);
end

str=sprintf('%0.3f',PSVi_median);
Header = strcat('Median PSVi_M_a_v_g = ',str);
title(Header);
MakePlotAxesNice


%%
% 6. PSVindex time series; iqr PSVindex

PSVi_iqr = iqr(PSVindex);

ax1 = nexttile;
hist(ax1,PSVindex,100)
xlim([0,1.1*max(PSVindex)]);

str=sprintf('%0.3f',PSVi_iqr);
Header = strcat('IQR PSVi_M_a_v_g_ = ',str);
title(Header);
xlabel('PSVi_M_a_v_g'); ylabel('N')
MakePlotAxesNice


%%
% Save figure

title_text = strcat(Title_names(Model)," Lat = ",str1," Long = ",str2);
title(t,title_text);
SaveFilePath = char(cell2mat(strcat(CurrentDirectory2,'\', PlotFileName)));
savefig(SaveFilePath)

