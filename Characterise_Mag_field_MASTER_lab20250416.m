% Script to plot magnetic field behaviour and characteristics for time
% series of gauss coefficients. Can be applied to outputs of geomagnetic
% field models, geodynamo sims, etc

%Required:

% -  a numeric matrix called gh_full of 120 rows and T columns where T is
% the number of timesteps. Each column is a single timestep and each row a
% gauss coefficient (in nT for observational models) starting at degree 1, order 0 in row 1 and
% degree 10, order 10 (g and h) in rows 119 and 120.

% - a numeric column vector called timesteps containing T evenly-speced time
%  values relating to columns in gh_full
% The units of this for observation models should be ka but should go
% forward in time from high to low values
% The units of this for dynamo sims should be magnetic diffusion times and
% should go forward in time 
%

% - The above saved as gh_full.mat and timesteps.mat in the working
% directory. This directory should be unique to that set of gauss
% coefficients as all outputs will also be saved there

% - CurrentDirectory specified as a the folder directly above the
% workingdirectory containing the file nm.mat

%- Set Model_List to 1, 2, or 3 (see below)

% - Set Model_names as a string (or string vector for multiple) containing
%  the name of the folder containing gh_full.mat and timesteps.mat
%
% - Set Title_names as a string (or string vector of same length as
% Model_names) containing the model name as you want it to appear on the
% plots


%Set Ref_SiteLat and Ref_SiteLat according to location of where you want
% your local time series to be l

clear;
delete(findall(0,'Type','figure'));
load('nm.mat') %gives degree and order numbers up to 10

%***Set constants
max_degreePSV = 10;  % Max degree and order for PSV analysis; default = 10
max_degreeTAF = 4;   % Max degree and order for TAF analysis; default = 4
Nt = 500; %Number of timesteps for all calcs of S and VDM_median etc; default = 500
Ref_SiteLat = 45; Ref_SiteLong = 0; %Reference Location for any time series of VGPs or VDMs etc.

%***Set Model Choice
Model_List = 1; %set set to 1 for single/multiple Observation-based models; 2 for single/multiple dynamo sims; 3 for list of paper dynamo sims

%****Set resolution
Maxlm_TAF = 10; %This only changes MAP. The parameters will be caluculated ot degree 4 regardless
%Always uses degree and order 10 for PSV

%****Set height for maps and TAF calcs
height = 0; %Earth's surface
%height = -2889; %CMB

%****Set central longitude for map
Long = -180:180; clong=1;
%Long(1,1:270) = [90:359]; Long(1,271:361) =[0:90]; clong=181;

%%

if Model_List == 1
    %load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Geomag_models_Datasets\Obs_Model_names.mat')
    
    gh_CorrFac = 1; %***CONVERT TO nT IF NECESSARY
    gh_CorrFac_DM = 1; gh_CorrFac_VDM = 1e-9; % Use this for OBS models to get VDM and Dipole Moment in correct units
    CurrentDirectory = "C:\Users\murray98\Documents\Bruce Buffett model\code\";
    %CurrentDirectory = "C:\Users\andyb\Dropbox\WORK_AB\PhD_Students and PDRAs\Mary Murray\Buffett model\";
    Model_names = [%"IGRF14"
        %"model017"; "model018"; "model019"; "model020";
        %"model021"; "model022"; "model023"; "model024";
        %"model025"; "model026"; 
        %"model027";
        %"model028"; "model029"; "model030";
        %"model031"; "model032"; "model033";
        %"model034"; "model035"; "model036";
        %"model037"; "model038";
        %"model104"; "model105"; "model106"; "model107"; "model108"; "model109"; 
        "model038a"; "model038";
        
        ]; %"Real1_10ka";%["Model_GGF100k";"Model_GGFSS70"];
    Title_names = [%"IGRF14"
        %"model017"; "model018"; "model019"; "model020";
        %"model021"; "model022"; "model023"; "model024";
        %"model025"; "model026"; 
        %"model027";
        %"model028"; "model029"; "model030";
        %"model031"; "model032"; "model033";
        %"model034"; "model035"; "model036";
        %"model037"; "model038";
        "model038a"; "model038";
        ];%["GGF100k";"GGFSS70"];
    %Model_names = Obs_Model_names(13,1);%[4 13 16 17])%(10:12,:)%(1:9,:)
    %Title_names = Obs_Model_names(13,2);
    Time_axis_label = "Time (ka)";

elseif Model_List == 2
    load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Geodynamo_Sims\Model_names.mat')
    gh_CorrFac = 1e6; %***CONVERT TO nT IF NECESSARY
    gh_CorrFac_DM = 1e6; gh_CorrFac_VDM = 1e-3; % Use this for DYNAMO models to get VDM and Dipole Moment in correct units
    CurrentDirectory = 'C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Geodynamo_Sims\';
    %CurrentDirectory = 'C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Papers\2018 06 PSVparameter\2019 08_submission1\Models\';
    
    Model_names = "LEDJ007"
    Title_names = "LEDJ007"
    %Model_names = Model_names(263)%:12)%:22);%(64,:)%%(62:67,:)%11:55,:)
    %Title_names = Model_names;
    Time_axis_label = "Time (kyr)";

elseif Model_List == 3
    %load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Papers\2020 12 RILs\20230715_Het_CMB_paper\Table_S1\TableS1_details.mat')
    %load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Papers\2020 12 RILs\20230715_Het_CMB_paper\ThermalHet.mat');
    %load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Papers\2020 12 RILs\20230715_Het_CMB_paper\ThermalHom.mat');
    %load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Papers\2020 12 RILs\20230715_Het_CMB_paper\TCHET1.mat');
    load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Papers\2020 12 RILs\20230715_Het_CMB_paper\THOM2_THET2_3.mat');
    %Model_names = Model_names(7,1);%(5:6,1);
    
    CurrentDirectory = 'C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Geodynamo_Sims\';
    gh_CorrFac = 1e6; %***CONVERT TO nT IF NECESSARY
    gh_CorrFac_DM = 1e6; gh_CorrFac_VDM = 1e-3; % Use this for DYNAMO models to get VDM and Dipole Moment in correct units
    Title_names = Model_names%(7,1);%(5:6,2);
    Time_axis_label = "Time (kyr)";
    
end
Nmodels = size(Model_names,1);
%%
for Model = 1: Nmodels

    clearvars -except n m Model_names Model_List Title_names CurrentDirectory Nmodels ...
        Model max_degreePSV max_degreeTAF Nt Ref_SiteLat Ref_SiteLong Time_axis_label gh_CorrFac_DM ...
        gh_CorrFac_VDM gh_CorrFac ResultsTable

    ModelName = Model_names(Model)
    CurrentDirectory2 = strcat(CurrentDirectory,ModelName);

    ghfullFilePath = strcat(CurrentDirectory2,'\gh_full.mat');
    ghfullFilePath = cell2mat(ghfullFilePath);
    load(char(ghfullFilePath))

    gh_full=gh_full .* gh_CorrFac; %***CONVERT TO nT IF NECESSARY

    timestepsfullFilePath = strcat(CurrentDirectory2,'\timesteps.mat');
    timestepsfullFilePath = cell2mat(timestepsfullFilePath);
    load(char(timestepsfullFilePath))

    ModelLength =size(timesteps,1);
    
    if Model_List == 1
        timesteps_dim = timesteps;
    else %geodynamo sim
        timesteps_dim = timesteps - min(timesteps);
        timesteps_dim = timesteps_dim ./ 5 .*1e3 ; %Converts from Magnetic diffusion time to kyr assuming 1 x MDF 200 kyr
    end
    
    ResultsTable.ModelDuration(Model,1) = (max(abs(timesteps_dim)) - min(abs(timesteps_dim)));

    %**********INSERT SCRIPTS OR FUNCTION HERE***************

    %AssessPSV %Generates Assess_PSV.mat FullResults.mat which gives S at even longitude and latitude spacing (non-uniform with density of S measurements increasing with latitude)
    %AssessPSV_VDM %Generates

     Characterise_Mag_field_global_time_series
     Characterise_Mag_field_local_time_series
     AssessPSV 
    %Power_Spectra
    %Characterise_Mag_field_ADF_NADF_time_series
    Characterise_Mag_field_ADF_NADF_No_plots

    ResultsTable.Model(Model,1) = ModelName;
    
    ResultsTable.Model_G_a_S_nocut(Model, 1) = ModelG_Results(1);
    ResultsTable.Model_G_b_S_nocut(Model, 1) = ModelG_Results(2);
    ResultsTable.Model_G_a_Svd(Model, 1) = ModelG_Results(4);
    ResultsTable.Model_G_b_Svd(Model, 1) = ModelG_Results(5);
    ResultsTable.Model_G_a_SpSE(Model, 1) = ModelG_Results(7);
    ResultsTable.Model_G_b_SpSE(Model, 1) = ModelG_Results(8);
    
    
    
    %%transistion time fraction (as defined by Sprain, 2019)
    ResultsTable.t_trans(Model, 1) = sum(d_colat >= 45 & d_colat <= 135) / numel(d_colat); 
    
    ResultsTable.VDM_median(Model,1) = median(VDM);
    ResultsTable.VDM_IQR(Model, 1) = iqr(VDM);

end

Results = table(ResultsTable.Model, ResultsTable.ModelDuration, ...
    ResultsTable.ADF_median, ResultsTable.ADF_iqr, ...
    ResultsTable.NADF_median, ResultsTable.NADF_iqr, ...
    ResultsTable.ADF_StabilisationTime, ResultsTable.NADF_StabilisationTime, ...
    ResultsTable.Model_G_a_S_nocut, ...
    ResultsTable.Model_G_b_S_nocut, ...
    ResultsTable.Model_G_a_Svd, ...
    ResultsTable.Model_G_b_Svd, ...
    ResultsTable.Model_G_a_SpSE, ...
    ResultsTable.Model_G_b_SpSE, ...
    ResultsTable.t_trans, ...
    ResultsTable.VDM_median, ResultsTable.VDM_IQR,...
    'VariableNames', ["Model", "Duration_kyr", "ADF_median", "ADF_iqr", ...
    "NADF_median", "NADF_iqr", "ADF_StableTime", ...
    "NADF_StableTime", ... 
    "Model_G_a_S_nocut", "Model_G_b_S_nocut", ...
    "Model_G_a_Svd", "Model_G_b_Svd", ...
    "Model_G_a_SpSE", "Model_G_b_SpSE", ...
    "t_trans", "VDM_median", "VDM_IQR" ,...
    ]);