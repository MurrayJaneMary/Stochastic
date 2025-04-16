% Plot:
% 1. VADM timeseries and median ADM
% 2. VADM distribution and iqr
% 3. VGP latitude timeseries; median VGP colat (polarity normalised);
% 4. VGP latitude distibution; iqr VGP colat (polarity normalised);
% 5. PSVindex time series; median PSVindex
% 6. PSVindex time series; iqr PSVindex


PlotFileName = "Global_TimeSeries_Plots_Lat.fig";

figure;
%set(gcf, 'Position',  [100, 100, 500, 250.*Nmodels]) %set final fig to 250 x Nmodels
fig = gcf;
fig.WindowState = 'maximized';
t=tiledlayout(3,2);

%%
% 1. Axial Dipole Moment timeseries and median ADM
Nt = size(timesteps_dim,1);

gh_full_ADM = gh_full .* gh_CorrFac_DM;
AxialDipoleMoment = DipoleMomentFunc(gh_full_ADM(1,:),0,0)';
ADM_median = median(AxialDipoleMoment);
ADM_median_line = [ADM_median ADM_median];

ax1= nexttile; %xlim([0,max(timesteps_dim)]);
plot(ax1,timesteps_dim,AxialDipoleMoment);hold on;
plot(ax1, [timesteps_dim(1,1), timesteps_dim(Nt,1)] , ADM_median_line);hold off;
xlim([min(timesteps_dim),max(timesteps_dim)]);
xlabel(Time_axis_label);

if timesteps_dim(1) > timesteps_dim(2)
    plot(ax1,flip(timesteps_dim),flip(AxialDipoleMoment));hold on;
    plot(ax1, flip([timesteps_dim(1,1), timesteps_dim(Nt,1)]) , ADM_median_line);hold off;
    set(gca, 'XDir', 'reverse'); xlabel('Time (ka)');
end

ylim([0 1.1*max(AxialDipoleMoment)]);
ylabel('Axial Dipole Moment (ZAm^2)');


str=sprintf('%0.2f',ADM_median);
Header = strcat('ADM_m_e_d_i_a_n = ',str);
title(Header);
MakePlotAxesNice

%%
% 2. DM distribution and iqr
ADM_iqr = iqr(AxialDipoleMoment);

ax1 = nexttile;
hist(ax1,AxialDipoleMoment,100)
xlim([0,1.1*max(AxialDipoleMoment)]);

str=sprintf('%0.2f',ADM_iqr);
Header = strcat('ADM_i_q_r = ',str);
title(Header);xlabel('Axial Dipole Moment (ZAm^2)'); ylabel('N')
MakePlotAxesNice
%%
% 3. Dipole tilt timeseries
g10 = gh_full(1,:)';
g11 = gh_full(2,:)';
h11 = gh_full(3,:)';
[dplat] = rad2deg(atan2(g10,sqrt(g11.^2 + h11.^2))) .* -1;

gh_full_n = polarity_norm_gauss(gh_full);
g10_n = gh_full_n(1,:)';
g11_n = gh_full_n(2,:)';
h11_n = gh_full_n(3,:)';
[d_colat] = 90-rad2deg(atan2(g10,sqrt(g11.^2 + h11.^2))) .* -1;
dipole_colat_median = median(d_colat);
dipole_lat_median_line = [90-dipole_colat_median 90-dipole_colat_median];

ax1= nexttile;
plot(ax1, timesteps_dim,dplat,'Color','k')
xlim([min(timesteps_dim),max(timesteps_dim)]); ylim([-90 90]);yticks(-90:45:90);
hold on
plot(ax1, [timesteps_dim(1,1), timesteps_dim(Nt,1)] , dipole_lat_median_line);
hold off
xlabel(Time_axis_label);

if timesteps_dim(1) > timesteps_dim(2)

    plot(ax1,flip(timesteps_dim),flip(dplat));hold on;
    plot(ax1, flip([timesteps_dim(1,1), timesteps_dim(Nt,1)]) , dipole_lat_median_line);hold off;
    set(gca, 'XDir', 'reverse'); xlabel('Time (ka)');
end

ylabel('pole latitude (°)');
%title('Dipole Tilt')
str=sprintf('%0.2f',dipole_colat_median);
Header = strcat('Dipole Colat_m_e_d_i_a_n = ',str);
title(Header);
MakePlotAxesNice
%%
% 4. Dipole Tilt distibution; iqr pole colat
Dipole_Colat_iqr = iqr(d_colat);

ax1 = nexttile;
hist(ax1,d_colat,100)

str=sprintf('%0.2f',ADM_iqr);
Header = strcat('Dipole Colat_i_q_r = ',str);
title(Header);xlabel('Dipole Colatitude (°)'); ylabel('N')
MakePlotAxesNice

%%
% 5. axial dipole and nonaxial dipole time series at CMB and time-average faxdip and fdip

%Find Axial dipole power at CMB
gh_full_AD = gh_full; %.* 1e-3; %change from nT to uT
gh_full_AD(2:3,:)=0;
[MagPower_AD_CMB, MeanField_AD_CMB] = magPowerCalc(gh_full_AD(1:3,:), 1, 1);

%Find total power at CMB
[MagPower_CMB, MeanField] = magPowerCalc(gh_full, 1, 10);

MaxRow = size(MagPower_CMB,1); MaxColumn = size(MagPower_CMB,2);

%MagField_CMB = sqrt(MagPower_CMB);

for timestep = 1: MaxColumn
    MagField_CMB_AxDipole(timestep) = sqrt(sumsqr(MeanField_AD_CMB(1,timestep)));
    MagField_CMB_Dipole(timestep) = sqrt(sumsqr(MeanField(1,timestep)));
    MagField_CMB_nondipole(timestep) = sqrt(sumsqr(MeanField(2:MaxRow,timestep)));
    MagField_CMB(timestep) = sqrt(sumsqr(MeanField(1:MaxRow,timestep)));
end

%MagPower_CMB_Dipole = sum(MagPower_CMB(1:3,:),1);
%MagPower_CMB_nondipole = sum(MagPower_CMB(4:MaxRow,:),1);
[fdip,fdip_timesteps] = fdip_calc(gh_full,max_degreePSV);
[faxdip,faxdip_timesteps] = faxdip_calc(gh_full,max_degreePSV);

ax1= nexttile;
plot(ax1, timesteps_dim,MagField_CMB_AxDipole,'Color','b','LineWidth',1); hold on;
plot(ax1, timesteps_dim,MagField_CMB_Dipole,'Color','k','LineWidth',1); hold on;
plot(ax1, timesteps_dim,MagField_CMB_nondipole,'Color','r','LineWidth',1);%hold off;
plot(ax1, timesteps_dim,MagField_CMB,'Color','g','LineWidth',0.5);hold off;

%plot(ax1, timesteps_dim,fdip_timesteps,'Color','r');
xlim([min(timesteps_dim),max(timesteps_dim)]); ylim([0,1.1*max(MagField_CMB)]);

if timesteps_dim(1) > timesteps_dim(2)
    plot(ax1, timesteps_dim,MagField_CMB_AxDipole,'Color','b','LineWidth',1); hold on;
    plot(ax1, flip(timesteps_dim), flip(MagField_CMB_Dipole),'Color','k','LineWidth',1);hold on;
    plot(ax1, flip(timesteps_dim), flip(MagField_CMB_nondipole),'Color','r','LineWidth',1);
    plot(ax1, flip(timesteps_dim), flip(MagField_CMB),'Color','g','LineWidth',0.5);hold off;
    set(gca, 'XDir', 'reverse'); xlabel('Time (ka)');
end

str1=sprintf('%0.2f',faxdip);
str2=sprintf('%0.2f',fdip);
Header = strcat('CMB: Time-averaged faxdip = ',str1, '; fdip = ',str2);
title(Header);xlabel(Time_axis_label); ylabel('field');
legend("Axial Dipole","Dipole","NonDipole","Total");
MakePlotAxesNice
%%
% 6. AD and NAD time series at surface and AD/NADmedian
gh_full_AD = gh_full .* 1e-3; %change from nT to uT
MagPower_AD = magPowerCalc(gh_full_AD(1:3,:), 0, 1);

gh_full_NAD = gh_full .* 1e-3; %change from nT to uT;
gh_full_NAD(1,:)=0;
MagPower_NAD = magPowerCalc(gh_full_NAD, 0, max_degreePSV);
MagPower_NAD = sum(MagPower_NAD,1);

for TimeStep = 1: Nt
    %TimeStep
    gh = gh_full(:,TimeStep);
    [~, ~, AD_NADSurf(1,TimeStep), ~, ~, ~]= EnergyFunc2(gh,n,m);
end

Median_AD_NADSurf = median(AD_NADSurf);
str=sprintf('%0.1f',Median_AD_NADSurf);

ax1= nexttile;
semilogy(ax1, timesteps_dim,MagPower_AD,'Color','k','LineWidth',1); hold on;
semilogy(ax1, timesteps_dim,MagPower_NAD,'Color','r');hold off;
xlim([min(timesteps_dim),max(timesteps_dim)]);
xlabel(Time_axis_label);

if timesteps_dim(1) > timesteps_dim(2)
    semilogy(ax1, flip(timesteps_dim), flip(MagPower_AD),'Color','k','LineWidth',1); hold on;
    semilogy(ax1, flip(timesteps_dim), flip(MagPower_NAD),'Color','r');hold off;
    set(gca, 'XDir', 'reverse'); xlabel('Time (ka)');
end

Header = strcat('Surface: Median AD/NAD = ',str);
title(Header); ylabel("Power");
legend("Axial Dipole","Non Axial Dipole");
MakePlotAxesNice

%%
% Save figure
title(t,Title_names(Model));
SaveFilePath = char(cell2mat(strcat(CurrentDirectory2,'\', PlotFileName)));
savefig(SaveFilePath)

