%Updated 20191024 to replace PSVn with ModelGfast_v2

% gh=nan(size(gh_full,1),1);
% AxDipSurf=nan(1,ModelLength); AxDip=nan(1,ModelLength); AD_NADSurf=nan(1,ModelLength);
% plat=nan(1,ModelLength); plong=nan(1,ModelLength);
% for TimeStep = 1: ModelLength
%     TimeStep
%     gh = gh_full(:,TimeStep);
% 
%     [AxDipSurf(1,TimeStep), AxDip(1,TimeStep), AD_NADSurf(1,TimeStep), ~, ~, ~]= EnergyFunc2(gh,n,m);
% 
%     [~,~, ~, ~, ~, ~, plat(1,TimeStep), ~, VDM(1,TimeStep), ~]...
%         =QueryModel(gh_full, TimeStep,slat,slong,max_degree);
% end

figure; set(gca, 'Fontsize', 24);
% subplot(1,3,1); 

% plot(timesteps_dim,gh_full(1,1:ModelLength)); title(Model_names(Model));
% xlabel(Time_axis_label); ylabel('g10 (arbitrary units)')
% 
% subplot(2,3,2); plot(timesteps_dim,plat); title('Site location = 45N 0E');
% xlabel(Time_axis_label); ylabel('VGP Latitude')
% platPerm=plat;
% 
% %Median_Odd_Even = median(Odd_Even);
% Median_AD_NADSurf = median(AD_NADSurf);
% 
% str=sprintf('%0.1f',Median_AD_NADSurf)
% Header = strcat('Median AD/NAD = ',str)
% subplot(2,3,3); semilogy(timesteps_dim,AD_NADSurf); title(Header);
% xlabel(Time_axis_label); ylabel('AD/NAD')
% 


%*********FINISHED WITH TIME SERIES*************
LongVector = [0:20:340]; %Set longitudes
LatVector =[-85:10:85]; %Set latitudes
Nd(1:Nt,1) = 999; kdir(1:Nt,1) = 99999;

%rand('state',sum(100*clock))

NVGPs = size(LongVector,2) .* size(LatVector,2) .* Nt
VGPs = nan(NVGPs,4); %VGPs = [dec, inc, VGPlat, VGPlong]

FullResults = nan(NVGPs./Nt, 9); 

for run=1:1
    
    NumLat = 1; NumLong = 1; location =1;
    
    row=1;
    while NumLong <= size(LongVector,2)
        slong = LongVector(NumLong)
        while NumLat <= size(LatVector,2)
            slat = LatVector(NumLat);
            
            
            [timestep, Sitelat, Sitelong, dec, inc, ~, plat, plong, ~, ~]...
                = RandDynamoPick(gh_full, Nt, 1, ModelLength, slat, slong);
            
            [platn, plongn, decn, incn] = FlipPolFunc(plat, plong, dec, inc);
            %Have VGPs and directions
            
            
            
            %Find S and S_vD
            [slatm, slongm, Decm, Incm, Maglatm, platm, plongm, Norig, S_nocut, Nnew, cutoff, low95, Svd, high95]...
                =VGPs_vD_func_NEW_vRKB(Sitelat, Sitelong, decn,incn, Nd, kdir, 0);
            PallatfromFishPole = 90-angdis([slat, slong, platm, plongm]);
            
            PropTrans = (Norig - Nnew)/Norig;
            
            
             [CentralLat, CentralLong, S_SE, Sprime_SE] = SphExpFunc2(platn, plongn); 
             PallatfromSEPole = 90-angdis([slat, slong, CentralLat, CentralLong]);
%             
             FullResults(location,:) = [location, slatm, slongm, PallatfromFishPole,...
                 S_nocut, PropTrans, Svd, PallatfromSEPole, Sprime_SE];
            
                        %test=test+1
            NumLat = NumLat + 1; location = location +1;
        end
        
        NumLong = NumLong +1; NumLat=1;
        
    end
    
end

%save('Assess_PSV.mat','FullResults')
PallatfromFishPole = FullResults(:,4); S_nocut = FullResults(:,5);
Svd = FullResults(:,7); PallatfromSEPole = FullResults(:,8); SpSE = FullResults(:,9);

[out_val, resid_nocut]...
    = ModelGfast_v2(PallatfromFishPole, S_nocut);
a_S_nocut=out_val(1,1); b_S_nocut=out_val(1,2);
PallatMOD2(1:181,1) = -90:1:90;
abestv2(1:181,1)=a_S_nocut; bbestv2(1:181,1)=b_S_nocut;
SbprimeMOD2_Snc = (abestv2.^2 + bbestv2.^2.*PallatMOD2.^2).^(1/2);


[out_val, resid_vd]...
    = ModelGfast_v2(PallatfromFishPole, Svd);
a_Svd=out_val(1,1); b_Svd=out_val(1,2);

PallatMOD2(1:181,1) = -90:1:90;
abestv2(1:181,1)=a_Svd; bbestv2(1:181,1)=b_Svd;
SbprimeMOD2_Svd = (abestv2.^2 + bbestv2.^2.*PallatMOD2.^2).^(1/2);

[out_val, resid_SE]...
    = ModelGfast_v2(PallatfromSEPole, SpSE);
a_SpSE=out_val(1,1); b_SpSE=out_val(1,2);

PallatMOD2(1:181,1) = -90:1:90;
abestv2(1:181,1)=a_SpSE; bbestv2(1:181,1)=b_SpSE;
SbprimeMOD2_SpSE = (abestv2.^2 + bbestv2.^2.*PallatMOD2.^2).^(1/2);

subplot(1,3,1); plot(FullResults(:,4),FullResults(:,5),'.'); hold on;
plot(PallatMOD2,SbprimeMOD2_Snc,'r-'); hold off
title('Fisher; no cutoff');
xlabel('Pallat'); ylabel('S'); xlim([-90 90]); ylim([0 80]);

subplot(1,3,2); plot(FullResults(:,4),FullResults(:,7),'.'); hold on;
plot(PallatMOD2,SbprimeMOD2_Svd,'r-'); hold off
title('Fisher; Vandamme');
xlabel('Pallat'); ylabel('S'); xlim([-90 90]); ylim([0 80]);

subplot(1,3,3); plot(FullResults(:,8),FullResults(:,9),'.'); hold on;
plot(PallatMOD2,SbprimeMOD2_SpSE,'r-'); hold off
title('Spherical Exponential');
xlabel('Pallat'); ylabel('Sprime'); xlim([-90 90]); ylim([0 80]);

%savefig('AssessPSV.fig')

ModelG_Results = [a_S_nocut,b_S_nocut,resid_nocut, a_Svd,b_Svd,resid_vd,...
    a_SpSE,b_SpSE,resid_SE];
ModelG_Results_real = real(ModelG_Results)

SaveFilePath = char(cell2mat(strcat(CurrentDirectory,'\Assess_PSV.mat')));
     save(SaveFilePath,'FullResults')
%    FullResults(location,:) = [location, slatm, slongm, PallatfromFishPole,...
%                 S_nocut, PropTrans, Svd, PallatfromSEPole, Sprime_SE]; 
     
    SaveFilePath = char(cell2mat(strcat(CurrentDirectory,'\Assess_PSV')));
    savefig(SaveFilePath)
%     
%     Summary(Model,1:4) = [size(gh_full,2), ModelLengthNDTU, StabilisationTimeNDTU, Nt]
    Summary(Model,:) = ModelG_Results_real;
    
     SaveFilePath = char(cell2mat(strcat(CurrentDirectory,'\ModelG_Results.mat')));
     save(SaveFilePath,'ModelG_Results_real')


%Make VDM time series plot

% max_yaxis = max(VDM) .* 1.1;
% figure
% scatter(timesteps_dim,VDM,'.')
% ylim([0 max_yaxis]);
% xlabel('Time (Myr)'); ylabel('VDM (nondimensionalised)')
