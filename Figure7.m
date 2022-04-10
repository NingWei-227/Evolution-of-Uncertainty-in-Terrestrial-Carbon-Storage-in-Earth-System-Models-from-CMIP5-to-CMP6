%% Fig4. Model-data comparison on global terrestrial NPP, GPP and CUE for two CMIPs
% Observational datasets
% GPP datasets:
% (1) MODIS17A2; 
% Unit: gC m-2 yr-1
% Running et al., 2015
% resolution£º 0.5 x 0.5 
%
% (2) GIMMSGPP
% reference: Smith NCC, 2016; available from: https://wkolby.org/data-code/
% unit: gC m-2 yr-1 
% resolution: 0.5 x 0.5
% standard estimation based on climate inputs from ECMWF, MERRA2 and NCEP
% 
% (3) FLUXCOM
% reference: Jung et al., 2017; https://www.bgc-jena.mpg.de/geodb/projects/Home.php
% Unit: gC m-2 day-1
% resolution: 0.5 x 0.5
%
% (4) VPM
% reference: Zhang et al., 2017;  https://doi.org/10.6084/m9.figshare.c.3789814
% Unit: gC m-2 yr-1
% resolution: 0.5 x 0.5

% NPP datasets used here are from MODIS17A2 and GIMMSGPP-NPP
% CUE (NPP/GPP) is calculated by using NPP and GPP data from the same source £¨MODIS17A2 and GIMMSGPP-NPP£©
clear;clc
%% load observational data
% GPP, Time: from 2001 to 2005;
% Data: FLUXCOM;  Unit: g C m-2 day-1
GPP_maps_1 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\FLUXCOM\yearlyGPP1980_2013.nc','GPP',[1 1 22],[Inf Inf 5]);
GPP_obs1(1:360,1:720,1:5) = NaN;
for i =1:5
    i
    GPPobs_M = GPP_maps_1(:,:,i);
    GPPobs_t = GPPobs_M';
    GPPobs_t = GPPobs_t.*365;  % convert unit to gC m-2 yr-1
    
    GPP_obs1(:,:,i) = GPPobs_t;   
end
% Data: MODIS17A2; Unit: gC m-2 yr-1
for i = 2001:2005
    [num,text,raw] = xlsread(['E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\MODIS17A2\MOD17A2GPPAnnual',num2str(i),'.csv']);
    for rowID=1:360
        for colID =1:720
            if raw{rowID,colID}=='NA'
                raw{rowID,colID}=nan;
            end  
        end
     end
    GPP_obs2(:,:,i-2000) = cell2mat(raw);
end
% Data: VPM; Unit: gC m-2 yr-1
for i = 2001:2005
    [num,text,raw] = xlsread(['E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\VPM\VPM_GPP_',num2str(i),'_Mean.csv']);
    for rowID=1:360
        for colID =1:720
            if raw{rowID,colID}=='NA'
                raw{rowID,colID}=nan;
            end  
        end
     end
    GPP_obs3(:,:,i-2000) = cell2mat(raw);
end
% Data: GIMMS-GPP
% reference: Smith NCC, 2016; available from: https://wkolby.org/data-code/
% unit: g C m-2 yr-1 
% resolution: 0.5 x 0.5
% standard estimation based on climate inputs from CRUNCEP-P1-standard-run
gpp_CRU = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\GIMMSNPP\GIMMSNPP\Annual\Annual\CRUNCEP\gpp_CRUNCEP_V4P1_Standard_1982_2016_Annual_GEO_30min.nc','GPP',[1 1 20],[Inf Inf 5]);
GPP_obs4(1:360,1:720,1:5) = NaN;
for i=1:5
    gppData4 = gpp_CRU(:,:,i)';
    GPP_obs4(:,:,i) = gppData4;
end

% merge GPP data
GPP_obs05_5yr(:,:,:,1) = GPP_obs1;
GPP_obs05_5yr(:,:,:,2) = GPP_obs2;
GPP_obs05_5yr(:,:,:,3) = GPP_obs3;
GPP_obs05_5yr(:,:,:,4) = GPP_obs4;

% unit: gC m-2 yr-1; 2001-2005 mean 
GPP_obs05(:,:,1) = nanmean(GPP_obs1,3);
GPP_obs05(:,:,2) = nanmean(GPP_obs2,3);
GPP_obs05(:,:,3) = nanmean(GPP_obs3,3);
GPP_obs05(:,:,4) = nanmean(GPP_obs4,3);

% load observational data of NPP
% Data: MODIS17A2;
% Unit: gC m-2 yr-1
% Running et al., 2015
% resolution£º 0.5 x 0.5
NPP_obs1(1:360,1:720,1:5) = NaN;
for i=2001:2005
    nppData = imread(['E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\MOD17A3_0.5_2016fool\MOD17A3_',num2str(i),'globalNPP.tif']);
    nppData = nppData.* 0.1; % Scale factor:0.1
    nppData(nppData<0) = NaN;
    nor_matrix(1:20,1:720) = NaN;
    sou_matrix(1:60,1:720) = NaN;
    
    nppMap = [nor_matrix; nppData; sou_matrix];
    NPP_obs1(:,:,i-2000) = nppMap;
end

% Data GIMMS-NPP
% reference: Smith NCC, 2016; available from: https://wkolby.org/data-code/
% unit: g C m-2 yr-1 
% resolution: 0.5 x 0.5
% standard estimation based on climate inputs from CRUNCEP-P1-standard-run
npp_CRU = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\GIMMSNPP\GIMMSNPP\Annual\Annual\CRUNCEP\npp_CRUNCEP_V4P1_Standard_1982_2016_Annual_GEO_30min.nc','NPP',[1 1 20],[Inf Inf 5]);
NPP_obs2(1:360,1:720,1:5) = NaN;
for i=1:5
    nppData2 = npp_CRU(:,:,i)';
    NPP_obs2(:,:,i) = nppData2;
end
% Data CARDAMOM estimates
% reference: Bloom et al., 2015 PNAS
% Unit: gC m-2 day-1
% resolution: 1 x 1
% decadal mean(2001-2010)
npp_CMOMnc = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\DS_10283_875\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_FL_NPP.nc','Mean');
npp_CMOM_re = npp_CMOMnc';
npp_CMOM1gre(1:180,1:360) = NaN;
for i=1:180
    npp_CMOM1gre(i,:) =  npp_CMOM_re(181-i,:);
end 
% regrided CARDAMOM data into 0.5x0.5 resolution
lat05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lat');
lon05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lon');
lat1 = 89.5:-1:-89.5; lat1 = lat1';
lon1 = -179.5:179.5;  lon1 = lon1';
[x05,y05] = meshgrid(lon05,lat05);
[x1,y1] = meshgrid(lon1,lat1);
npp_CMOM05 = interp2(x1,y1,npp_CMOM1gre,x05,y05,'linear');
% convert unit from gC m-2 day-1 into gC m-2 yr-1
NPP_obs3 = npp_CMOM05.*365;

% merge NPP data
% unit: gC m-2 yr-1; 
NPP_obs05_5yr(:,:,:,1) = NPP_obs1;
NPP_obs05_5yr(:,:,:,2) = NPP_obs2;

% unit: gC m-2 yr-1; 2001-2005 mean 
NPP_obs05(:,:,1) = nanmean(NPP_obs1,3);
NPP_obs05(:,:,2) = nanmean(NPP_obs2,3);
NPP_obs05(:,:,3) = NPP_obs3;

NPP_obs05KG = NPP_obs05.*10^(-3); % change unit from gC m-2 yr-1 into KgC m-2 yr-1
GPP_obs05KG = GPP_obs05.*10^(-3); % change unit from gC m-2 yr-1 into KgC m-2 yr-1

% observation-based estimates on CUE
% reference: Bloom et al., 2015 PNAS
% Unit: gC m-2 day-1
% resolution: 1 x 1
npp_CMOMnc = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\DS_10283_875\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_FL_NPP.nc','Mean');
npp_CMOM_25th = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\DS_10283_875\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_FL_NPP.nc','25th_percentile');
npp_CMOM_75th = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\DS_10283_875\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_FL_NPP.nc','75th_percentile');

gpp_CMOMnc = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\DS_10283_875\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_FL_GPP.nc','Mean');
gpp_CMOM_25th = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\DS_10283_875\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_FL_GPP.nc','25th_percentile');
gpp_CMOM_75th = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\DS_10283_875\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_FL_GPP.nc','75th_percentile');

npp_CMOM_re = npp_CMOMnc';
npp_CMOM25_re = npp_CMOM_25th';
npp_CMOM75_re = npp_CMOM_75th';

gpp_CMOM_re = gpp_CMOMnc';
gpp_CMOM25_re = gpp_CMOM_25th';
gpp_CMOM75_re = gpp_CMOM_75th';

npp_CMOM1gre(1:180,1:360) = NaN;
npp_CMOM1gre_25(1:180,1:360) = NaN;
npp_CMOM1gre_75(1:180,1:360) = NaN;
gpp_CMOM1gre(1:180,1:360) = NaN;
gpp_CMOM1gre_25(1:180,1:360) = NaN;
gpp_CMOM1gre_75(1:180,1:360) = NaN;
for i=1:180
    npp_CMOM1gre(i,:) =  npp_CMOM_re(181-i,:);
    npp_CMOM1gre_25(i,:) =  npp_CMOM25_re(181-i,:);
    npp_CMOM1gre_75(i,:) =  npp_CMOM75_re(181-i,:);
    
    gpp_CMOM1gre(i,:) =  gpp_CMOM_re(181-i,:);
    gpp_CMOM1gre_25(i,:) =  gpp_CMOM25_re(181-i,:);
    gpp_CMOM1gre_75(i,:) =  gpp_CMOM75_re(181-i,:);
end 
% regrided CARDAMOM data into 0.5x0.5 resolution
lat05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lat');
lon05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lon');
lat1 = 89.5:-1:-89.5; lat1 = lat1';
lon1 = -179.5:179.5;  lon1 = lon1';
[x05,y05] = meshgrid(lon05,lat05);
[x1,y1] = meshgrid(lon1,lat1);
npp_CMOM05 = interp2(x1,y1,npp_CMOM1gre,x05,y05,'linear');
npp_CMOM05_25 = interp2(x1,y1,npp_CMOM1gre_25,x05,y05,'linear');
npp_CMOM05_75 = interp2(x1,y1,npp_CMOM1gre_75,x05,y05,'linear');

gpp_CMOM05 = interp2(x1,y1,gpp_CMOM1gre,x05,y05,'linear');
gpp_CMOM05_25 = interp2(x1,y1,gpp_CMOM1gre_25,x05,y05,'linear');
gpp_CMOM05_75 = interp2(x1,y1,gpp_CMOM1gre_75,x05,y05,'linear');
% the spatial data of CUE
CUE_obs05_avg = npp_CMOM05./gpp_CMOM05;
CUE_obs05_25 = npp_CMOM05_25./gpp_CMOM05_25;
CUE_obs05_75 = npp_CMOM05_75./gpp_CMOM05_75;
CUE_obs05(:,:,1) = CUE_obs05_25;
CUE_obs05(:,:,2) = CUE_obs05_avg;
CUE_obs05(:,:,3) = CUE_obs05_75;
CUE_obs05(CUE_obs05>1) = NaN; % omit unreasonable values
CUE_obs05(CUE_obs05<0) = NaN;


% global estimation on CUE
area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv');
npp_CMOM05_pg = npp_CMOM05.*365.* area05.*10^(-15);    % Unit: PgC yr-1
npp_CMOM05_pgGB = nansum(npp_CMOM05_pg(:))
gpp_CMOM05_pg = gpp_CMOM05.*365.* area05.*10^(-15);    % Unit: PgC yr-1
gpp_CMOM05_pgGB = nansum(gpp_CMOM05_pg(:))

npp_CMOM05_25pg = npp_CMOM05_25.*365.* area05.*10^(-15);    % Unit: PgC yr-1
npp_CMOM05_25pgGB = nansum(npp_CMOM05_25pg(:))
gpp_CMOM05_25pg = gpp_CMOM05_25.*365.* area05.*10^(-15);    % Unit: PgC yr-1
gpp_CMOM05_25pgGB = nansum(gpp_CMOM05_25pg(:))

npp_CMOM05_75pg = npp_CMOM05_75.*365.* area05.*10^(-15);    % Unit: PgC yr-1
npp_CMOM05_75pgGB = nansum(npp_CMOM05_75pg(:))
gpp_CMOM05_75pg = gpp_CMOM05_75.*365.* area05.*10^(-15);    % Unit: PgC yr-1
gpp_CMOM05_75pgGB = nansum(gpp_CMOM05_75pg(:))

CUE_obsGB(1) = npp_CMOM05_pgGB./gpp_CMOM05_pgGB
CUE_obsGB(2) = npp_CMOM05_25pgGB./gpp_CMOM05_25pgGB
CUE_obsGB(3) = npp_CMOM05_75pgGB./gpp_CMOM05_75pgGB

clearvars -except NPP_obs05KG GPP_obs05KG CUE_obs05 NPP_obs05_5yr GPP_obs05_5yr CUE_obsGB

%% Global estimates on GPP, NPP and CUE
% CMIP5
cd('G:\CMIP5\4_MatData\temporal_data\hist_rcp85')
load('NPP_cmip5_tmp.mat')
load('GPP_cmip5_tmp.mat')
load('CUE_cmip5_tmp.mat')
% For HadGEM2-ES,GFDL-ESM2G and MRI-ESM1, their historical simulations did not fully cover the period of 1850-2005.
% We thus use NaN to supplement years without data  
NaNgf = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
NPPgf_tmp = [NaNgf; NPPgf_tmp]; 
GPPgf_tmp = [NaNgf; GPPgf_tmp];
CUEgf_tmp = [NaNgf; CUEgf_tmp];

NaNhad = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
NPPhad_tmp = [NaNhad; NPPhad_tmp];
GPPhad_tmp = [NaNhad; GPPhad_tmp];
CUEhad_tmp = [NaNhad; CUEhad_tmp];

NPPmri_tmp = [NaN; NPPmri_tmp];
GPPmri_tmp = [NaN; GPPmri_tmp];
CUEmri_tmp = [NaN; CUEmri_tmp];

NPP5_all(1:5,1:11) = NaN;
NPP5_all(:,1) = NPPbcc_tmp(152:156);
NPP5_all(:,2) = NPPcan_tmp(152:156);
NPP5_all(:,3) = NPPccsm_tmp(152:156);
NPP5_all(:,4) = NPPhad_tmp(152:156);
NPP5_all(:,5) = NPPipsl_tmp(152:156);
NPP5_all(:,6) = NPPmiroc_tmp(152:156);
NPP5_all(:,7) = NPPmpi_tmp(152:156);
NPP5_all(:,8) = NPPnor_tmp(152:156);
NPP5_all(:,9) = NPPbnu_tmp(152:156);
NPP5_all(:,10) = NPPgf_tmp(152:156);
NPP5_all(:,11) = NPPmri_tmp(152:156);

GPP5_all(1:5,1:11) = NaN;
GPP5_all(:,1) = GPPbcc_tmp(152:156);
GPP5_all(:,2) = GPPcan_tmp(152:156);
GPP5_all(:,3) = GPPccsm_tmp(152:156);
GPP5_all(:,4) = GPPhad_tmp(152:156);
GPP5_all(:,5) = GPPipsl_tmp(152:156);
GPP5_all(:,6) = GPPmiroc_tmp(152:156);
GPP5_all(:,7) = GPPmpi_tmp(152:156);
GPP5_all(:,8) = GPPnor_tmp(152:156);
GPP5_all(:,9) = GPPbnu_tmp(152:156);
GPP5_all(:,10) = GPPgf_tmp(152:156);
GPP5_all(:,11) = GPPmri_tmp(152:156);

% CMIP6
cd('H:\CMIP56_Csink\4_MatData\temporal_data\hist_ssp585')
load('2NPP_cmip6_tmp.mat')
load('2GPP_cmip6_tmp.mat')
load('2CUE_cmip6_tmp.mat')  
leg6_str = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'}    
NPP6_all = [NPPbcc_tmp(152:156),NPPcan_tmp(152:156),NPPcesm_tmp(152:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   NPPuk_tmp(152:156),NPPipsl_tmp(152:156),NPPmic_tmp(152:156),...            % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   NPPmpi_tmp(152:156),NPPnor_tmp(152:156),...                              % MPI-ESM1-2-LR, NorESM2
   NPPass_tmp(152:156),NPPcnrm_tmp(152:156),NPPec_tmp(152:156)];              % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg
GPP6_all = [GPPbcc_tmp(152:156),GPPcan_tmp(152:156),GPPcesm_tmp(152:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   GPPuk_tmp(152:156),GPPipsl_tmp(152:156),GPPmic_tmp(152:156),...          % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   GPPmpi_tmp(152:156),GPPnor_tmp(152:156),...                           % MPI-ESM1-2-LR, NorESM2
   GPPass_tmp(152:156),GPPcnrm_tmp(152:156),GPPec_tmp(152:156)];            % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

clearvars -except NPP_obs05KG GPP_obs05KG CUE_obs05 NPP_obs05_5yr GPP_obs05_5yr...
                  NPP5_sp05 NPP6_sp05 GPP5_sp05 GPP6_sp05 CUE5_sp05 CUE6_sp05 ...
                  NPP5_all GPP5_all NPP6_all GPP6_all CUE_obsGB
area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv');


NPP_obs05_Pg = NPP_obs05_5yr.* area05.*10^(-15);    % convert Unit into PgC yr-1
NPP_obs5yr_Pg = nansum(NPP_obs05_Pg,1);  NPP_obs5yr_Pg = nansum(NPP_obs5yr_Pg,2); 
NPP_obs5yr_Pg = squeeze(NPP_obs5yr_Pg)

GPP_obs05_Pg = GPP_obs05_5yr.* area05.*10^(-15);    % Unit: PgC yr-1
GPP_obs5yr_Pg = nansum(GPP_obs05_Pg,1);  GPP_obs5yr_Pg = nansum(GPP_obs5yr_Pg,2); 
GPP_obs5yr_Pg = squeeze(GPP_obs5yr_Pg)

%% Figure
% open figure window and set position
figure
set(gcf,'position',[100 100 940,531.2])
% Panel (b) and (c): Model-data comparison on CUE
% observational range of CUE
CUE_min = min(CUE_obsGB(:))
CUE_max = max(CUE_obsGB(:))

% CUE in CMIP5 and CMIP6
CUE5_5yr_range = NPP5_all./GPP5_all
CUE6_5yr_range = NPP6_all./GPP6_all
% 2001-2005 mean
CUE5_5yr_avg = nanmean(CUE5_5yr_range,1)
CUE6_5yr_avg = nanmean(CUE6_5yr_range,1)
CUE5_5yr_min = min(CUE5_5yr_range)
CUE6_5yr_min = min(CUE6_5yr_range)
CUE5_5yr_max = max(CUE5_5yr_range)
CUE6_5yr_max = max(CUE6_5yr_range)

% preparing data for plotting
ID_models = 1:11;
Markers = {'d','o','^','s','s','*','o','v','o','*','d'};
mycolor5 = [255 0 0; 153 51 255; 237 176 33; ...%BCC CanESM2 CCSM4
            0 197 205; 0 205 0; 207 194 124;...   %Had IPSL Miroc
            255 99 71; 65 105 255;...  %mpi NorM
            0 0 0; 158 131 149; 219 136 153]./255;  %BNU, GFDL, mri 
mycolor6 = [255 0 0; 153 51 255; 237 176 33;...    %BCC-CSM2-MR, CanESM5, CESM2
           0 197 205; 0 205 0; 207 194 124;...     %UKESM1-0-LL, IPSL-CM6A-LR,MIROC-ES2L
           255 99 71; 65 105 255;...               %MPI-ESM1-2-LR,NorESM-LM
            0 0 0; 158 131 149; 119 136 153]./255; %ACCESS-ESM1-5, CNRM-ESM2-1  EC-Earth3-Veg  


Models5 = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'};
Models6 = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'};
    
CUE5_plot = [ID_models', CUE5_5yr_avg', CUE5_5yr_min', CUE5_5yr_max']
CUE6_plot = [ID_models', CUE6_5yr_avg', CUE6_5yr_min', CUE6_5yr_max']
% sort data based on CUE 
CUE5_plot_sort = sortrows(CUE5_plot,2,'descend')   
CUE6_plot_sort = sortrows(CUE6_plot,2,'descend')

% Panel(b)
Panel_bc = tight_subplot(2,1,[0.17 0],[0.18 0.08],[0.68 0.02])
axes(Panel_bc(1))
hold on
X = 1:11
for i =1:11
   h5(i) = plot(X(i),CUE5_plot_sort(i,2),'Marker',Markers{i},...
       'MarkerEdgeColor', mycolor5(CUE5_plot_sort(i,1),:),...
       'MarkerSize',10,'LineStyle','none','LineWidth',2)    
end
set(gca, 'YLim',[0.3 0.6],'XLim',[0 12]);
x=[0 12 12 0];
y=[CUE_min CUE_min CUE_max CUE_max]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.3)
set(gca,'linewidth',1.2,'box','on')
set(gca,'Fontname','Arial','FontSize',10);
Ylabel = gca
set(Ylabel.YAxis,'FontSize',12)
set(gca,'XTickLabelRotation',40);
set(gca,'YTickLabelMode','auto');
ylabel('CUE','Fontname','Arial','FontSize',13)
xticks(1:11);
xticklabels({Models5{CUE5_plot_sort(1:11,1)}});
text(0.5,0.57,'(b) CMIP5','FontName','Arial','FontSize',12)

% Panel(c)
axes(Panel_bc(2))
hold on
X = 1:11
for i =1:11
   h6(i) = plot(X(i),CUE6_plot_sort(i,2),'Marker',Markers{i},...
       'MarkerEdgeColor', mycolor6(CUE6_plot_sort(i,1),:),...
       'MarkerSize',10,'LineStyle','none','LineWidth',2);   
end
set(gca, 'YLim',[0.3 0.6],'XLim',[0 12]);
x=[0 12 12 0];
y=[CUE_min CUE_min CUE_max CUE_max]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.3)
set(gca,'linewidth',1.2,'box','on')
set(gca,'Fontname','Arial','FontSize',10);
Ylabel = gca
set(Ylabel.YAxis,'FontSize',12)
set(gca,'XTickLabelRotation',40);
set(gca,'YTickLabelMode','auto');
ylabel('CUE','Fontname','Arial','FontSize',13)
xticks(1:11);
xticklabels({Models6{CUE6_plot_sort(1:11,1)}});
text(0.5,0.57,'(c) CMIP6','FontName','Arial','FontSize',12)


%% Panel(a): model-data comparison on GPP and NPP
Panel_a = tight_subplot(1,1,[0 0],[0.1 0.3],[0.08 0.52])
hold on
for i=1:11
    plot(GPP5_all(:,i),NPP5_all(:,i),'Marker',Markers{find(CUE5_plot_sort(:,1) == i)},...
       'MarkerEdgeColor', mycolor5(i,:),...
       'MarkerSize',10,'LineStyle','none','LineWidth',2)  
    plot(GPP6_all(:,i),NPP6_all(:,i),'Marker',Markers{find(CUE6_plot_sort(:,1) == i)},...
       'MarkerEdgeColor', mycolor6(i,:),...
       'MarkerSize',10,'LineStyle','none','LineWidth',2)  
    
end


set(gca, 'YLim',[40 100],'XLim',[80 260]);
% add the observational range
% Note that, GPP and NPP data from GIMMSGPP-NPP was calculated as ensemble mean 
GPP_global_Pg = nanmean(GPP_obs5yr_Pg,1)
% GPP range from the min to the max
GPP_obsMax = max(GPP_global_Pg)
GPP_obsMin = min(GPP_global_Pg)
% NPP
NPP_global_Pg = nanmean(NPP_obs5yr_Pg,1)
ITO_NPP = [50.6 68.4] % the ranges of observation-based estimates on NPP [mean+(-)SD] estimated from Ito et al.,2011, Table2, NPP over 2000s
% NPP range from the min to the max
NPP_obsMax = max([NPP_global_Pg ITO_NPP])
NPP_obsMin = min([NPP_global_Pg ITO_NPP])

% shading area for GPP 
y=[40 100 100 40];
x=[GPP_obsMin GPP_obsMin GPP_obsMax GPP_obsMax]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.3)
% shading area for NPP
x=[80 260 260 80];
y=[NPP_obsMin NPP_obsMin NPP_obsMax NPP_obsMax]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.6,0.6,0.6],'EdgeColor','none','FaceAlpha',0.3)
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
set(gca,'linewidth',1.2,'box','on')
set(gca,'Fontname','Arial','FontSize',12);
ylabel('NPP (PgC yr^-^1)','Fontname','Arial','FontSize',13)
xlabel('GPP (PgC yr^-^1)','Fontname','Arial','FontSize',13)
set(gca,'XLim',[80 260],'YLim',[40 100]);
%leg56 = legend({'CMIP5','CMIP6','GPP Obs','NPP Obs'})
text(92,95,'(a)','FontName','Arial','FontSize',12)

% GPP density plot for simulations from CMIP5 and CMIP6
Panel_Agpp = tight_subplot(1,1,[0 0],[0.705 0.06],[0.08 0.52])
hold on
[fGPP_cmip5, xiGPP_cmip5] = ksdensity(GPP5_all(:));
[fGPP_cmip6, xiGPP_cmip6] = ksdensity(GPP6_all(:));
plot(xiGPP_cmip5,fGPP_cmip5,'LineWidth',1.8,'color',[0.30,0.75,0.93]);
plot(xiGPP_cmip6,fGPP_cmip6,'LineWidth',1.8,'color',[1.00,0.07,0.65]);
set(gca,'XLim',[80 260]);
axis off
% NPP density polt for simulations from CMIP5 and CMIP6
Panel_Anpp = tight_subplot(1,1,[0 0],[0.1 0.3],[0.4807 0.38])
hold on
[fNPP_cmip5, xiNPP_cmip5] = ksdensity(NPP5_all(:));
[fNPP_cmip6, xiNPP_cmip6] = ksdensity(NPP6_all(:));
cmips(1) = plot(fNPP_cmip5,xiNPP_cmip5,'LineWidth',1.8,'color',[0.30,0.75,0.93]);
cmips(2) = plot(fNPP_cmip6,xiNPP_cmip6,'LineWidth',1.8,'color',[1.00,0.07,0.65]);
set(gca,'YLim',[38 102]);
axis off

leg_panel = legend(cmips,{'CMIP5','CMIP6'})
set(leg_panel,'FontName','Arial','FontSize',10,'Position',[0.355,0.607,0.098,0.066],...
    'color','w','EdgeColor','k')


