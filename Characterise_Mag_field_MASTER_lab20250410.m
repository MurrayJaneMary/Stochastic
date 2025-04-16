% Script to plot magnetic field behaviour and characteristics for time
% series of gauss coefficients. Can be applied to outputs of geomagnetic
% field models, geodynamo sims, etc

%Required:

% -  a numeric matrix called gh_full of 120 rows and T columns where T is
% the number of timesteps. Each column is a single timestep and each row a
% gauss coefficient starting at degree 1, order 0 in row 1 and
% degree 10, order 10 (g and h) in rows 119 and 120.

% - a numeric column vector called timesteps containing T evenly-speced time
%  values relating to columns in gh_full
% The units of this for observation models should be years
% The units of this for dynamo sims should be magnetic diffusion times
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
Model_List = 1; %set set to 1 for single/multiple Observation-based/stochastic models; 2 for single/multiple dynamo sims; 3 for list of paper dynamo sims

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
    %**************DEFINE CURRENT DIRECTORY AND FILE***********************
    CurrentDirectory = "C:\Users\murray98\Documents\Bruce Buffett model\code\"
    %CurrentDirectory = "C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Geomag_models_datasets\"
    %CurrentDirectory = "C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Geodynamo_Sims\"
    %Model_names = "Model_GGF100k"
    %Title_names = "Model GGF100k"
    Model_names = "model1"
    Title_names = "model1"
    %*************************************************************

    %load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Geomag_models_Datasets\Obs_Model_names.mat')
    gh_CorrFac_DM = 1; gh_CorrFac_VDM = 1e-6; % Use this for OBS models to get VDM and Dipole Moment in correct units

        %Model_names = Obs_Model_names(13,1);%[4 13 16 17])%(10:12,:)%(1:9,:)
    %Title_names = Obs_Model_names(13,2);
    Time_axis_label = "Time (yr)"
    


elseif Model_List == 2
    %**************DEFINE CURRENT DIRECTORY AND FILE***********************
    CurrentDirectory = "C:\Users\sligh\Dropbox\Analyse_SH_model_ouputs\LEDS023\"
    Model_names = "LEDS023"
    Title_names = "LEDS023"
    %*************************************************************

    gh_CorrFac_DM = 1e6; gh_CorrFac_VDM = 1e-3; % Use this for DYNAMO models to get VDM and Dipole Moment in correct units
    
    Time_axis_label = "Time (Myr)"
    




elseif Model_List == 3
    %load a specific list (string vector) of model names each giving the folder where the GCs are stored

    %**************DEFINE CURRENT DIRECTORY AND FILE***********************
    %CurrentDirectory = "3C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Geodynamo_Sims\"
    %load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Papers\2020 12 RILs\20230715_Het_CMB_paper\TableS1_details.mat')
    %load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Papers\2020 12 RILs\20230715_Het_CMB_paper\ThermalHet.mat');
    %load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Papers\2020 12 RILs\20230715_Het_CMB_paper\ThermalHom.mat');
    %load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Papers\2020 12 RILs\20230715_Het_CMB_paper\TCHET1.mat');
    %load('C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Papers\2020 12 RILs\20230715_Het_CMB_paper\THOM2_THET2_3.mat');
    %*************************************************************

        %Model_names = Model_names(7,1);%(5:6,1);

    CurrentDirectory = 'C:\Users\andyb\Dropbox\WORK_AB\Andy_Research\Geodynamo_Sims\';
    gh_CorrFac_DM = 1e6; gh_CorrFac_VDM = 1e-3; % Use this for DYNAMO models to get VDM and Dipole Moment in correct units
    Title_names = Model_names%
    Time_axis_label = "Time (Myr)"
    

end

Nmodels = size(Model_names,1);
%%
for Model = 1: Nmodels

    clearvars -except n m Model_names Model_List Title_names CurrentDirectory Nmodels ...
        Model max_degreePSV max_degreeTAF Nt Ref_SiteLat Ref_SiteLong Time_axis_label gh_CorrFac_DM ...
        gh_CorrFac_VDM Power_axis_label

    ModelName = Model_names(Model)
    CurrentDirectory2 = strcat(CurrentDirectory,ModelName);

    ghfullFilePath = strcat(CurrentDirectory2,'\gh_full.mat');
    ghfullFilePath = cell2mat(ghfullFilePath);
    load(char(ghfullFilePath))

    timestepsfullFilePath = strcat(CurrentDirectory2,'\timesteps.mat');
    timestepsfullFilePath = cell2mat(timestepsfullFilePath);
    load(char(timestepsfullFilePath))

    ModelLength =size(timesteps,1);

    if Model_List == 1
        timesteps_dim = timesteps;
    else %geodynamo sim
        timesteps_dim = timesteps - min(timesteps);
        timesteps_dim = timesteps_dim ./ 5; %Converts from Magnetic diffusion time to Myr assuming 1 x MDF 200 kyr
    end

    %**********CHOOSE SCRIPTS HERE***************
    
    Characterise_Mag_field_global_time_series
    Characterise_Mag_field_local_time_series
    AssessPSV 
    Power_Spectra

end