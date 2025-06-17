% Plot:
% 1. time series of rms Axial_Dipole_Field and rms Non-Axial_Dipole_Field at Earth's surface 
% 2. Histogram of rms Axial_Dipole_Field and Median and IQR
% 3. Histogram of rms non-Axial_Dipole_Field and Median and IQR
% and averaging times

PlotFileName = "ADF_NADF_TimeSeries_Plots.fig";

figure;
%set(gcf, 'Position',  [100, 100, 500, 250.*Nmodels]) %set final fig to 250 x Nmodels
fig = gcf;
fig.WindowState = 'maximized';
t=tiledlayout(3,2);

Threshold = 0.01; %Fraction of leniancy for time-averaged median 
%%

% 1. axial dipole and nonaxial dipole time series at surface and
% time-average ADF/NADF

%Find Axial dipole power at CMB
[ADF,NADF] = NADF_calc(gh_full,max_degreePSV);
ADF = ADF ./ 1000; %convert to uT
NADF = NADF ./ 1000; %convert to uT

%Calculate running medians
[ADF_RunningMedian N_ADF N_ADF_Sufficient] =...
    RunningMedianThreshold(timesteps_dim, ADF, Threshold);
[NADF_RunningMedian N_NADF N_NADF_Sufficient] =...
    RunningMedianThreshold(timesteps_dim, NADF, Threshold);


%Make a sampling of timesteps to improve speed

if size(timesteps_dim,1) > 10000
    step = round(size(timesteps_dim,1)./100);
    timesteps_dim_sample(1:100,1)= timesteps_dim(1:step:100*step);
    gh_full_sample = gh_full(1:15,1:step:100*step);
else
    timesteps_dim_sample = timesteps_dim;
    gh_full_sample = gh_full;
end

[term_RunningMedian(1:15,:) N_terms(1:15,:) N_terms_Sufficient(1:15,:)] =...
    RunningMedianThreshold(timesteps_dim_sample, gh_full_sample(1:15,:), Threshold);

%Normalise individual terms by their max values
max_gh = max(abs(gh_full),[],2);
term_RunningMedianMAXNORM = term_RunningMedian(1:15,:) ./ max_gh(1:15);



%Find time after which median is always within threshold

ADF_StabilisationTime = abs(timesteps_dim(1) - timesteps_dim(max(find(N_ADF_Sufficient==0))+1)); 
if isempty(ADF_StabilisationTime) == 1  
    ADF_StabilisationTime = 0;
end
NADF_StabilisationTime = abs(timesteps_dim(1) - timesteps_dim(max(find(N_NADF_Sufficient==0))+1)); 
if isempty(NADF_StabilisationTime) == 1  
    NADF_StabilisationTime = 0;
end

for term = 1:15
term_StabilisationTime(term,1) = abs(timesteps_dim(1) - timesteps_dim(max(find(N_terms_Sufficient(term,:)==0))+1)); 
end



MaxRMS = max(max([ADF; NADF]));

ax1= nexttile([1 2]);
plot(ax1, timesteps_dim,ADF,'Color','b','LineWidth',1); hold on;
plot(ax1, timesteps_dim,NADF,'Color','r','LineWidth',1); 
plot(ax1, timesteps_dim,ADF_RunningMedian,'Color','k','LineWidth',0.7); 
plot(ax1, timesteps_dim,NADF_RunningMedian, '--','Color','k','LineWidth',0.7);
hold off;

%plot(ax1, timesteps_dim,fdip_timesteps,'Color','r');
xlim([min(timesteps_dim),max(timesteps_dim)]); ylim([0,1.1*MaxRMS]);

if timesteps_dim(1) > timesteps_dim(2)
     set(gca, 'XDir', 'reverse'); 
end

%ADF_StabilisationTime = ADF_StabilisationTimeAbs/max(abs(timesteps_dim));
%NADF_StabilisationTime = NADF_StabilisationTimeAbs/max(abs(timesteps_dim));

xlabel(Time_axis_label);
ylabel('RMS Field strength (uT)');



% 
% str1=sprintf('%0.2f',faxdip);
% str2=sprintf('%0.2f',fdip);
% Header = strcat('CMB: Time-averaged faxdip = ',str1, '; fdip = ',str2);
%title(Title_names(Model));
legend("Axial Dipole","Nonaxial-Dipole");
MakePlotAxesNice

%%
% 1a. Time series of all terms with l and m < 4 normalised by their max
% values

ax1= nexttile([1 2]);
plot(ax1, timesteps_dim_sample,term_RunningMedianMAXNORM(1,:),'Color','b','LineWidth',2); hold on;
for i=2:15
    plot(ax1, timesteps_dim_sample,term_RunningMedianMAXNORM(i,:),'LineWidth',0.5);
end

if timesteps_dim(1) > timesteps_dim(2)
     set(gca, 'XDir', 'reverse'); 
end

yline(ax1, 0, '--', 'LineWidth', 1.5); % Add dashed line at y = 0
xlim([min(timesteps_dim),max(timesteps_dim_sample)]);
legendLabels = ["g10", "g11","h11","g20","g21","h21","g22","h22","g30","g31","h31","g32","h32","g33","h33"];
legend(ax1, legendLabels);
xlabel(Time_axis_label);
ylabel('Median / Max field strengths');
MakePlotAxesNice

%%
% 2. Histogram of rms Axial_Dipole_Field and Median and IQR
ADF_median = median(ADF);
ADF_iqr = iqr(ADF);

ax1 = nexttile;
hist(ax1,ADF,100)
xlim([0,1.1*MaxRMS]);

str1=sprintf('%0.1f',ADF_median);
str2=sprintf('%0.1f',ADF_iqr);
Header = strcat('ADF_m_e_d_i_a_n = ', str1, ' uT; ADM_i_q_r = ',str2, ' uT');
title(Header);xlabel('RMS Axial Dipole Field (uT)'); ylabel('N')
MakePlotAxesNice

%%
% 3. Histogram of rms NonAxial_Dipole_Field and Median and IQR
NADF_median = median(NADF);
NADF_iqr = iqr(NADF);

ax1 = nexttile;
hist(ax1,NADF,100)
xlim([0,1.1*MaxRMS]);

str1=sprintf('%0.1f',NADF_median);
str2=sprintf('%0.1f',NADF_iqr);
Header = strcat('NADF_m_e_d_i_a_n = ', str1, ' uT; NADM_i_q_r = ',str2, ' uT');
title(Header);xlabel('RMS Nonaxial Dipole Field (uT)'); ylabel('N')
MakePlotAxesNice
%%
% Save figure and data
title(t,Title_names(Model));
SaveFilePath = char(cell2mat(strcat(CurrentDirectory2,'\', PlotFileName)));
savefig(SaveFilePath)

ResultsTable.ADF_median(Model,1) = ADF_median;
ResultsTable.ADF_iqr(Model,1) = ADF_iqr;
ResultsTable.NADF_median(Model,1) = NADF_median;
ResultsTable.NADF_iqr(Model,1) = NADF_iqr;
ResultsTable.ADF_StabilisationTime(Model,1) = ADF_StabilisationTime;
ResultsTable.NADF_StabilisationTime(Model,1) = NADF_StabilisationTime;