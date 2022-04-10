clear
clc

%% CMIP5
% veg and  soil carbon(accounting for cLitter and cSoil)
%  CMIP5 cVeg data at 1x1 resolution 
cd('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP5')
cVeg5_BCC =  ncread('1gre_cPool3_yr_bcc-csm1-1-m_r1i1p1_185001-210012.nc','cVeg');
cVeg5_CAN = ncread('1gre_cPool3_yr_CanESM2_r1i1p1_185001-210012.nc','cVeg');
cVeg5_CCSM = ncread('1gre_cVeg_yr_CCSM4_r1i1p1_185001-210012.nc','cVeg');
cVeg5_HAD = ncread('1gre_cPool2_yr_HadGEM2-ES_r1i1p1_186001-210012.nc','cVeg');

cVeg5_IPSL = ncread('1gre_cPool3_yr_IPSL-CM5A-MR_r1i1p1_185001-210012.nc','cVeg');
cVeg5_MIROC = ncread('1gre_cPool3_yr_MIROC-ESM_r1i1p1_185001-210012.nc','cVeg');
cVeg5_MPI = ncread('1gre_cPool3_yr_MPI-ESM-MR_r1i1p1_185001-210012.nc','cVeg');
cVeg5_NOR = ncread('3_1gre_yr_cPool3_mon_NorESM1-M_r1i1p1_185001-210012.nc','cVeg');

cVeg5_BNU = ncread('1gre_cPool3_yr_BNU-ESM_r1i1p1_185001-210012.nc','cVeg');
cVeg5_GF = ncread('1gre_cPool2_yr_GFDL-ESM2G_r1i1p1_186101-210012.nc','cVeg');
cVeg5_MRI = ncread('1gre_cPool3_yr_MRI-ESM1_r1i1p1_185101-210012.nc','cVeg');

% CMIP5 clitter and SOM data at 1x1 resolution 
cLit5_BCC =  ncread('1gre_cPool3_yr_bcc-csm1-1-m_r1i1p1_185001-210012.nc','cLitter');
cLit5_CAN = ncread('1gre_cPool3_yr_CanESM2_r1i1p1_185001-210012.nc','cLitter');
cLit5_CCSM = ncread('1gre_cLitter_yr_CCSM4_r1i1p1_185001-210012.nc','cLitter');
%cLit_HAD = ncread('1gre_cPool2_yr_HadGEM2-ES_r1i1p1_186001-210012.nc','cLitter');
cLit5_IPSL = ncread('1gre_cPool3_yr_IPSL-CM5A-MR_r1i1p1_185001-210012.nc','cLitter');
cLit5_MIROC = ncread('1gre_cPool3_yr_MIROC-ESM_r1i1p1_185001-210012.nc','cLitter');
cLit5_MPI = ncread('1gre_cPool3_yr_MPI-ESM-MR_r1i1p1_185001-210012.nc','cLitter');
cLit5_NOR = ncread('3_1gre_yr_cPool3_mon_NorESM1-M_r1i1p1_185001-210012.nc','cLitter');
cLit5_BNU = ncread('1gre_cPool3_yr_BNU-ESM_r1i1p1_185001-210012.nc','cLitter');
%cLit_GF = ncread('1gre_cPool2_yr_GFDL-ESM2G_r1i1p1_186101-210012.nc','cLitter');
cLit5_MRI = ncread('1gre_cPool3_yr_MRI-ESM1_r1i1p1_185101-210012.nc','cLitter');

cSOM5_BCC =  ncread('1gre_cPool3_yr_bcc-csm1-1-m_r1i1p1_185001-210012.nc','cSoil');
cSOM5_CAN = ncread('1gre_cPool3_yr_CanESM2_r1i1p1_185001-210012.nc','cSoil');
cSOM5_CCSM = ncread('1gre_cSoil_yr_CCSM4_r1i1p1_185001-210012.nc','cSoil');
cSOM5_HAD = ncread('1gre_cPool2_yr_HadGEM2-ES_r1i1p1_186001-210012.nc','cSoil');
cSOM5_IPSL = ncread('1gre_cPool3_yr_IPSL-CM5A-MR_r1i1p1_185001-210012.nc','cSoil');
cSOM5_MIROC = ncread('1gre_cPool3_yr_MIROC-ESM_r1i1p1_185001-210012.nc','cSoil');
cSOM5_MPI = ncread('1gre_cPool3_yr_MPI-ESM-MR_r1i1p1_185001-210012.nc','cSoil');
cSOM5_NOR = ncread('3_1gre_yr_cPool3_mon_NorESM1-M_r1i1p1_185001-210012.nc','cSoil');
cSOM5_BNU = ncread('1gre_cPool3_yr_BNU-ESM_r1i1p1_185001-210012.nc','cSoil');
cSOM5_GF = ncread('1gre_cPool2_yr_GFDL-ESM2G_r1i1p1_186101-210012.nc','cSoil');
cSOM5_MRI = ncread('1gre_cPool3_yr_MRI-ESM1_r1i1p1_185101-210012.nc','cSoil');


cSoil5_BCC=cLit5_BCC+cSOM5_BCC;   cSoil5_CAN = cLit5_CAN+cSOM5_CAN;        cSoil5_CCSM = cLit5_CCSM + cSOM5_CCSM;
cSoil5_HAD=cSOM5_HAD;             cSoil5_IPSL = cLit5_IPSL + cSOM5_IPSL;   cSoil5_MIROC = cLit5_MIROC + cSOM5_MIROC;
cSoil5_MPI=cLit5_MPI+cSOM5_MPI;   cSoil5_NOR = cLit5_NOR + cSOM5_NOR;      cSoil5_BNU = cLit5_BNU + cSOM5_BNU;
cSoil5_GF=cSOM5_GF;               cSoil5_MRI = cLit5_MRI + cSOM5_MRI;
              
HAD_NaN(1:360,1:180,1:10) = NaN;
cSoil5_HAD = cat(3,HAD_NaN,cSoil5_HAD);
cVeg5_HAD = cat(3,HAD_NaN,cVeg5_HAD);

GF_NaN(1:360,1:180,1:11) = NaN;
cSoil5_GF = cat(3,GF_NaN,cSoil5_GF);
cVeg5_GF = cat(3,GF_NaN,cVeg5_GF);

MRI_NaN(1:360,1:180,1) = NaN;
cSoil5_MRI = cat(3,MRI_NaN,cSoil5_MRI);
cVeg5_MRI = cat(3,MRI_NaN,cVeg5_MRI);     

%% CMIP6
% CMIP6 cVeg data              
cd('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6')
cVeg6_BCC = ncread('2_1gre_cVeg_yr_Lmon_BCC-CSM2-MR_historical_r1i1p1f1_gn_1850-2014.nc','cVeg');
cVeg6_CAN = ncread('3_1gre_cVeg_yr_CanESM5_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');
cVeg6_CESM = ncread('4_1gre_cVeg_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');
cVeg6_UK = ncread('11_1gre_cVeg_yr_UKESM1-0-LL_historical_r1i1p1f2_gn_185001-201412.nc','cVeg');
cVeg6_IPSL = ncread('7_1gre_cPool3_yr_IPSL-CM6A-LR_historical_r1i1p1f1_gr_185001-201412.nc','cVeg');
cVeg6_MIC = ncread('8_1gre_cVeg_yr_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','cVeg');
cVeg6_MPI = ncread('9_1gre_cVeg_yr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');
cVeg6_NOR = ncread('10_1gre_cVeg_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');
cVeg6_ASS = ncread('1_1gre_yr_cVeg_Lmon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');
cVeg6_CNRM = ncread('5_1gre_yr_cVeg_Lmon_CNRM-ESM2-1_historical_r1i1p1f2_gr_185001-201412.nc','cVeg');
cVeg6_EC = ncread('6_1gre_cVeg_yr_EC-Earth3-Veg_historical_r1i1p1f1_gr_185001-201412.nc','cVeg');

% CMIP6 clitter and SOM data at 1x1 resolution 
cLit6_BCC = ncread('2_1gre_cLitter_yr_Lmon_BCC-CSM2-MR_historical_r1i1p1f1_gn_1850-2014.nc','cLitter');
cLit6_CAN = ncread('3_1gre_cLitter_yr_CanESM5_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');
cLit6_CESM = ncread('4_1gre_cLitter_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');
cLit6_IPSL = ncread('7_1gre_cPool3_yr_IPSL-CM6A-LR_historical_r1i1p1f1_gr_185001-201412.nc','cLitter');
cLit6_MIC = ncread('8_1gre_cLitter_yr_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','cLitter');
cLit6_MPI = ncread('9_1gre_cLitter_yr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');
cLit6_NOR = ncread('10_1gre_cLitter_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');
cLit6_ASS = ncread('1_1gre_yr_cLitter_Lmon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');
cLit6_CNRM = ncread('5_1gre_yr_cLitter_Lmon_CNRM-ESM2-1_historical_r1i1p1f2_gr_185001-201412.nc','cLitter');
cLit6_EC = ncread('6_1gre_cLitter_yr_EC-Earth3-Veg_historical_r1i1p1f1_gr_185001-201412.nc','cLitter');

% CESM2 and NorESM2 simulated soil carbon storage along soil profiles
% For benchmark analysis, we used cSoilAbove1m for the two models.
% evaluation of the full-depth SOM was conducted elsewhere
cSOM6_BCC = ncread('2_1gre_cSoil_yr_Emon_BCC-CSM2-MR_historical_r1i1p1f1_gn_1850-2014.nc','cSoil');
cSOM6_CAN = ncread('3_1gre_cSoil_yr_CanESM5_historical_r1i1p1f1_gn_185001-201412.nc','cSoil');
cSOM6_CESM_1m = ncread('4_1gre_cSoilAbove1m_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cSoilAbove1m');
cSOM6_UK = ncread('11_1gre_cSoil_yr_UKESM1-0-LL_historical_r1i1p1f2_gn_185001-201412.nc','cSoil');
cSOM6_IPSL = ncread('7_1gre_cPool3_yr_IPSL-CM6A-LR_historical_r1i1p1f1_gr_185001-201412.nc','cSoil');
cSOM6_MIC = ncread('8_1gre_cSoil_yr_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','cSoil');
cSOM6_MPI = ncread('9_1gre_cSoil_yr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_185001-201412.nc','cSoil');
cSOM6_NOR_1m = ncread('10_1gre_cSoilAbove1m_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cSoilAbove1m');
cSOM6_ASS = ncread('1_1gre_yr_cSoil_Emon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc','cSoil');
cSOM6_CNRM = ncread('5_1gre_yr_cSoil_Emon_CNRM-ESM2-1_historical_r1i1p1f2_gr_185001-201412.nc','cSoil');
cSOM6_EC = ncread('6_1gre_cSoil_yr_EC-Earth3-Veg_historical_r1i1p1f1_gr_185001-201412.nc','cSoil');

cSoil6_BCC = cLit6_BCC+cSOM6_BCC;    cSoil6_CAN = cLit6_CAN+cSOM6_CAN;      cSoil6_CESM = cLit6_CESM+cSOM6_CESM_1m;
cSoil6_UK = cSOM6_UK;                cSoil6_IPSL = cLit6_IPSL+cSOM6_IPSL;   cSoil6_MIC = cLit6_MIC+cSOM6_MIC;
cSoil6_MPI = cLit6_MPI+cSOM6_MPI;    cSoil6_NOR = cLit6_NOR+cSOM6_NOR_1m;      cSoil6_ASS = cLit6_ASS+cSOM6_ASS;
cSoil6_CNRM = cLit6_CNRM+cSOM6_CNRM; cSoil6_EC = cLit6_EC+cSOM6_EC;

%%
% cLand (cVeg and cSoil) from 2001 to 2005 
% unit:KgC m-2
% CMIP5
sp5_cLand_kgC(:,:,:,1) = cVeg5_BCC(:,:,152:156) + cSoil5_BCC(:,:,152:156);
sp5_cLand_kgC(:,:,:,2) = cVeg5_CAN(:,:,152:156) + cSoil5_CAN(:,:,152:156);
sp5_cLand_kgC(:,:,:,3) = cVeg5_CCSM(:,:,152:156) + cSoil5_CCSM(:,:,152:156);
sp5_cLand_kgC(:,:,:,4) = cVeg5_HAD(:,:,152:156) + cSoil5_HAD(:,:,152:156);
sp5_cLand_kgC(:,:,:,5) = cVeg5_IPSL(:,:,152:156) + cSoil5_IPSL(:,:,152:156);
sp5_cLand_kgC(:,:,:,6) = cVeg5_MIROC(:,:,152:156) + cSoil5_MIROC(:,:,152:156);
sp5_cLand_kgC(:,:,:,7) = cVeg5_MPI(:,:,152:156) + cSoil5_MPI(:,:,152:156);
sp5_cLand_kgC(:,:,:,8) = cVeg5_NOR(:,:,152:156) + cSoil5_NOR(:,:,152:156);
sp5_cLand_kgC(:,:,:,9) = cVeg5_BNU(:,:,152:156) + cSoil5_BNU(:,:,152:156);
sp5_cLand_kgC(:,:,:,10) = cVeg5_GF(:,:,152:156) + cSoil5_GF(:,:,152:156);
sp5_cLand_kgC(:,:,:,11) = cVeg5_MRI(:,:,152:156)+ cSoil5_MRI(:,:,152:156);
% CMIP6
sp6_cLand_kgC(:,:,:,1) = cVeg6_BCC(:,:,152:156) + cSoil6_BCC(:,:,152:156);
sp6_cLand_kgC(:,:,:,2) = cVeg6_CAN(:,:,152:156) + cSoil6_CAN(:,:,152:156);
sp6_cLand_kgC(:,:,:,3) = cVeg6_CESM(:,:,152:156) + cSoil6_CESM(:,:,152:156);
sp6_cLand_kgC(:,:,:,4) = cVeg6_UK(:,:,152:156) + cSoil6_UK(:,:,152:156);
sp6_cLand_kgC(:,:,:,5) = cVeg6_IPSL(:,:,152:156) + cSoil6_IPSL(:,:,152:156);
sp6_cLand_kgC(:,:,:,6) = cVeg6_MIC(:,:,152:156) + cSoil6_MIC(:,:,152:156);
sp6_cLand_kgC(:,:,:,7) = cVeg6_MPI(:,:,152:156) + cSoil6_MPI(:,:,152:156);
sp6_cLand_kgC(:,:,:,8) = cVeg6_NOR(:,:,152:156) + cSoil6_NOR(:,:,152:156);
sp6_cLand_kgC(:,:,:,9) = cVeg6_ASS(:,:,152:156) + cSoil6_ASS(:,:,152:156);
sp6_cLand_kgC(:,:,:,10) = cVeg6_CNRM(:,:,152:156) + cSoil6_CNRM(:,:,152:156);
sp6_cLand_kgC(:,:,:,11) = cVeg6_EC(:,:,152:156) + cSoil6_EC(:,:,152:156);

% 2001-2005 mean cLand
cLand11_map5_kgC = nanmean(sp5_cLand_kgC,3); cLand11_map5_kgC = squeeze(cLand11_map5_kgC);      % unit: KgC m-2
cLand11_map6_kgC = nanmean(sp6_cLand_kgC,3); cLand11_map6_kgC = squeeze(cLand11_map6_kgC);      % unit: KgC m-2

% convert to global maps
cLand5_map11_kgC(1:180,1:360,1:11) = NaN;
cLand6_map11_kgC(1:180,1:360,1:11) = NaN;

for i =1:11
    i
    
    cLand5_avg = cLand11_map5_kgC(:,:,i);
    cLand6_avg = cLand11_map6_kgC(:,:,i);
    
    cLand5_t = cLand5_avg';
    cLand6_t = cLand6_avg';
    
    cLand5_LR(1:180,1:360) = NaN;
    cLand6_LR(1:180,1:360) = NaN;
    
    cLand5_LR(:,1:180) = cLand5_t(:,181:360);
    cLand5_LR(:,181:360) = cLand5_t(:,1:180);
    cLand6_LR(:,1:180) = cLand6_t(:,181:360);
    cLand6_LR(:,181:360) = cLand6_t(:,1:180);
    
    map_cLand5(1:180,1:360) = NaN;
    map_cLand6(1:180,1:360) = NaN;
    for k=1:180
        map_cLand5(181-k,:) = cLand5_LR(k,:);
        map_cLand6(181-k,:) = cLand6_LR(k,:);
    end
    
    cLand5_map11_kgC(:,:,i) = map_cLand5;
    cLand6_map11_kgC(:,:,i) = map_cLand6;      
end

% convert to 0.5x0.5 degree
cLand5_map11_05kgC(1:360,1:720,1:11) = NaN;
cLand6_map11_05kgC(1:360,1:720,1:11) = NaN;

lat05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lat');
lon05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lon');
lat1 = 89.5:-1:-89.5; lat1 = lat1';
lon1 = -179.5:179.5;  lon1 = lon1';

[x05,y05] = meshgrid(lon05,lat05);
[x1,y1] = meshgrid(lon1,lat1);

for models=1:11
    cLand5_M = interp2(x1,y1,cLand5_map11_kgC(:,:,models),x05,y05,'linear');
    cLand6_M = interp2(x1,y1,cLand6_map11_kgC(:,:,models),x05,y05,'linear');
    
    cLand5_map11_05kgC(:,:,models) = cLand5_M;
    cLand6_map11_05kgC(:,:,models) = cLand6_M;
end

mask5 = cLand5_map11_05kgC(:,:,1);
mask5(~isnan(mask5)) = 1;
mask6 = cLand6_map11_05kgC(:,:,1);
mask6(~isnan(mask6)) = 1;

cLand5_map11_05kgC = cLand5_map11_05kgC.*mask5;
cLand6_map11_05kgC = cLand6_map11_05kgC.*mask6;

%% Global terrestrial C storage in CMIP5 and CMIP6
load E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\4_MatData\cVeg_cSoil_CMIP56_native.mat
cVeg5_org_avg11 = nanmean(cVeg5_org_PG);
cSoil5_org_avg11 = nanmean(cSoil5_org_PG);
cVeg6_org_avg11 = nanmean(cVeg6_org_PG);
cSoil6_org_avg11 = nanmean(cSoil6_org_PG);

cLand5_Pg_11M = cVeg5_org_avg11 + cSoil5_org_avg11;
cLand6_Pg_11M = cVeg6_org_avg11 + cSoil6_org_avg11; 

clearvars -except  cLand5_map11_05kgC cLand6_map11_05kgC cLand5_Pg_11M cLand6_Pg_11M mask5 mask6

%% load observational data
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\Biomass\ORNL\Global_Maps_C_Density_2010_1763\cVeg05_KgCm2.mat')
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\HWSD\HWSD_1247\HWSDv2_05.mat')
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\LandGIS\LandGIS05_1m.mat')
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\SoilGrids\SoilGrid05_1mkg.mat')
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\Circumpolar_cSoil\NCSCDv2_Circumpolar_netCDF_05deg\NCSCD05_1mkg.mat')
  
cSoil_data(:,:,1) = HWSD1m_05;
cSoil_data(:,:,2) = LandGIS_1mKg;
cSoil_data(:,:,3) = SoilGrid_1mkg;
cSoil_data(:,:,4) = NCSCDgb_1mKg;

cSoil_mean = nanmean(cSoil_data,3);
cLand_obs_mean = cVeg_obs05kg + cSoil_mean;
cLand_obs_range = cSoil_data + cVeg_obs05kg;

cLand_obs_min(1:360,1:720) = NaN;
cLand_obs_max(1:360,1:720) = NaN;
for i=1:360
    for j= 1:720            
        cLand_obs_min(i,j) = nanmin(cLand_obs_range(i,j,:));
        cLand_obs_max(i,j) = nanmax(cLand_obs_range(i,j,:));          
    end
end

polar_mask = NCSCDgb_1mKg;
polar_mask(~isnan(polar_mask)) = 1;

% omit deserts and other regions where GPP<0.01 kgC m-2 yr-1 based on the  criteria of Fan et al., 2019
% Observational dataset of tuaE 
tuaE_map_obs = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','tau');
tuaE_median1m_obs = tuaE_map_obs(:,:,1)';   % considering 0-1m soil C
tuaE_median1m_obs(~isnan(tuaE_median1m_obs)) = 1;
Mask_tua = tuaE_median1m_obs;

cSoil_mean = cSoil_mean.* Mask_tua;
cSoil_data = cSoil_data.* Mask_tua;

cVeg_obs05kg = cVeg_obs05kg.* Mask_tua;
cVeg_obs_un05kg = cVeg_obs_un05kg.* Mask_tua;

cLand_obs_mean = cLand_obs_mean.* Mask_tua;
cLand_obs_range = cLand_obs_range.* Mask_tua;
cLand_obs_min  = cLand_obs_min.* Mask_tua;
cLand_obs_max = cLand_obs_max.* Mask_tua;

cLand5_map11_05kgC = cLand5_map11_05kgC.*Mask_tua;
cLand6_map11_05kgC = cLand6_map11_05kgC.*Mask_tua;

cLand5_map11_05kgC = cLand5_map11_05kgC.*mask5;
cLand6_map11_05kgC = cLand6_map11_05kgC.*mask6;


clearvars -except cSoil_mean cSoil_data ...
                  cVeg_obs05kg cVeg_obs_un05kg ...
                  cLand_obs_mean cLand_obs_range cLand_obs_min  cLand_obs_max...
                  polar_mask Mask_tua cLand5_map11_05kgC cLand6_map11_05kgC cLand5_Pg_11M cLand6_Pg_11M mask5 mask6

%% Figure
Models5 = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'}              
Models6 = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'} 
mycolor5 = [255 0 0; 153 51 255; 237 176 33; ...%BCC CanESM2 CCSM4
            0 197 205; 0 205 0; 207 194 124;...   %Had IPSL Miroc
            255 99 71; 65 105 255;...  %mpi NorM
            0 0 0; 158 131 149; 119 136 153]./255;  %BNU, GFDL, mri 
mycolor6 = [255 0 0; 153 51 255; 237 176 33;...    %BCC-CSM2-MR, CanESM5, CESM2
           0 197 205; 0 205 0; 207 194 124;...     %UKESM1-0-LL, IPSL-CM6A-LR,MIROC-ES2L
           255 99 71; 65 105 255;...               %MPI-ESM1-2-LR,NorESM-LM
            0 0 0; 158 131 149; 119 136 153]./255; %ACCESS-ESM1-5, CNRM-ESM2-1  EC-Earth3-Veg      
    
    
X5_endT = table('Size',[11 3],'VariableType',{'string','double','double'},...
    'VariableName',{'Moldes','X','Nlimitation'});
X6_endT = table('Size',[11 3],'VariableType',{'string','double','double'},...
    'VariableName',{'Moldes','X','Nlimitation'});

X5_endT.Moldes = Models5';
X6_endT.Moldes = Models6';
X5_endT.X = cLand5_Pg_11M';
X6_endT.X = cLand6_Pg_11M';

NL5 = []; NL6 = []; % number of models considered nutrient limitation
NL5 = [0 0 1 0 0 0 0 1 0 0 0]';
NL6 = [0 0 1 1 1 1 1 1  1 0 1]';
X5_endT.Nlimitation = NL5;
X6_endT.Nlimitation = NL6;

avgX5 = nanmean(X5_endT.X); avgX6 = nanmean(X6_endT.X);
sdX5 = nanstd(X5_endT.X); sdX6 = nanstd(X6_endT.X);

%%
cLand5_SD11_05kgC = nanstd(cLand5_map11_05kgC,0,3); 
cLand6_SD11_05kgC = nanstd(cLand6_map11_05kgC,0,3);

cLand5_SD11_05kgC(cLand5_SD11_05kgC>999) = NaN;
cLand6_SD11_05kgC(cLand6_SD11_05kgC>999) = NaN;

load H:\CMIP56\Figure_Codes2\Codes_case1\spatial_Codes\MyColor16.mat
load H:\CMIP56\CMIP6_trace\CMIP6_Models\1_Mat_data\5_Models_spatial\MyColor10_X2.mat
load H:\CMIP56\CMIP6_trace\CMIP6_Models\1_Mat_data\5_Models_spatial\MyColor10_SD2.mat
load('E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure2_GPP\mycolor_GPPsd.mat');
%ax = gca;
%mycolor_sd = colormap(ax)
%save ('E:\1_Mycase\3_CMIP56_Cland\5_Version4\Figure1\mycolor_sd.mat','mycolor_sd') 
load 'E:\1_Mycase\3_CMIP56_Cland\5_Version4\Figure1\mycolor_sd.mat'

%CMIP5
cLand5_SD_M11 = cLand5_SD11_05kgC;
cLand5_SD_M11(302:360,:) = [];
cLand5_SD_M11 = flipud(cLand5_SD_M11);
raster5_cLandSD = georasterref('RasterSize',size(cLand5_SD_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);


figure
set(gcf,'position',[100 100 870,540])
maps = tight_subplot(2,1,[-0.08 -0.085],[0.25 -0.018],[-0.01 0.49])
fig_cmip5 = maps(1)
axes(fig_cmip5)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
setm(fig_cmip5, 'FFaceColor', [0.68,0.88,0.96])
framem('FLineWidth',1,'FEdgeColor','none')
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3,'DefaultEdgeColor',[0.4 0.4 0.4])
framem('FLineWidth',1)
h = geoshow(cLand5_SD_M11,raster5_cLandSD, 'DisplayType','surface','Zdata',zeros(size(cLand5_SD_M11)),'CData',cLand5_SD_M11);
colormap(fig_cmip5,mycolor_sd)
caxis([0 25])
set(gca,'box','off')
axis off
colorbar('off')
text(-2.312,2.0883, '(a) CMIP5','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);

% CMIP6
cLand6_SD_M11 = cLand6_SD11_05kgC;
cLand6_SD_M11(302:360,:) = [];
cLand6_SD_M11 = flipud(cLand6_SD_M11);
raster6_cLandSD = georasterref('RasterSize',size(cLand6_SD_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);
fig_cmip6 = maps(2)
axes(fig_cmip6)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
setm(fig_cmip6, 'FFaceColor', [0.68,0.88,0.96])
framem('FLineWidth',1,'FEdgeColor','none')
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3,'DefaultEdgeColor',[0.4 0.4 0.4])
framem('FLineWidth',1)
h = geoshow(cLand6_SD_M11,raster6_cLandSD, 'DisplayType','surface','Zdata',zeros(size(cLand6_SD_M11)),'CData',cLand6_SD_M11);
colormap(fig_cmip6,mycolor_sd)
caxis([0 25])
set(gca,'box','off')
axis off
colorbar('off')
text(-2.312,2.0883, '(d) CMIP6','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
h_bar5 = colorbar
h_bar5.Location = 'southoutside'
h_bar5.Position = [0.085,0.28,0.32714351425943,0.0298];
h_bar5.FontName = 'Arial'
h_bar5.FontSize = 10;
text(0.0183,-2.3264,'Standard deviation in cLand',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11)

obs_avg = 2285; % Table 3
obs_min = 1516; % Table 3
obs_max = 2835; % Table 3
pane1 = tight_subplot(2,1,[0.15 0.25],[0.32 0.15],[0.082 0.85]) 
axes(pane1(1))
hold on
for i =1:11
    CMIP5_lines(i) = plot([0.5,1],[X5_endT.X(i) X5_endT.X(i)] ,'LineWidth',3,'LineStyle','-','color',mycolor5(i,:))
    
end
plot(1.3, obs_avg,'*','MarkerSize',8,'color','k')
errorbar(1.3, obs_avg, obs_avg-obs_min, obs_max-obs_avg,'color','k','LineWidth',1)
text(1.27,1320,'Obs','FontName','Arial','FontSize',10,'color','r','Rotation',-5)

%H_pa5 = patch([0.75,0.75,1.25,1.25],[avgX5-sdX5,avgX5+sdX5,avgX5+sdX5,avgX5-sdX5],[0.6,0.6,0.6]);
%set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
%avg5 = line([0.65 1.25],[avgX5 avgX5],'color','k','linewidth',2);
set(gca,'XLim',[0.4 1.6],'YLim',[500 3600],'LineWidth',1.2,'box','off')
set(gca,'Xcolor','none')
set(gca,'Fontname','Arial','FontSize',10)
set(gca,'YTickLabelMode','auto');
ylabel('Pg C','Fontname','Arial','FontSize',11)
set(gca,'color','none','TickDir','out','TickLength',[0.03 0.035])

leg5 = legend(CMIP5_lines,Models5,'Color','none',...
    'FontName','Arial','FontSize',8,...
    'NumColumns',3,'Position',[0.0797,0.038,0.4299,0.142],...
    'EdgeColor',[0.3 0.3 0.3],'LineWidth',1)
title(leg5,'CMIP5 models','FontName','Arial','FontSize',10)




% CMIP6
axes(pane1(2)) 
hold on
for i =1:11
    CMIP6_lines(i) = plot([0.5,1],[X6_endT.X(i) X6_endT.X(i)] ,'LineWidth',3,'LineStyle','-','color',mycolor6(i,:))
    
end
plot(1.3, obs_avg,'*','MarkerSize',8,'color','k')
errorbar(1.3, obs_avg, obs_avg-obs_min, obs_max-obs_avg,'color','k','LineWidth',1)
text(1.27,1320,'Obs','FontName','Arial','FontSize',10,'color','r','Rotation',-5)

set(gca,'XLim',[0.4 1.6],'YLim',[500 3600],'LineWidth',1.2,'box','off')
set(gca,'Xcolor','none')
set(gca,'Fontname','Arial','FontSize',10)
set(gca,'YTickLabelMode','auto');
ylabel('Pg C','Fontname','Arial','FontSize',11)
set(gca,'color','none','TickDir','out','TickLength',[0.03 0.035])

leg6 = legend(CMIP6_lines,Models6,'Color','none',...
    'FontName','Arial','FontSize',8,...
    'NumColumns',3,'Position',[0.5363,0.038,0.4557,0.142],...
    'EdgeColor',[0.3 0.3 0.3],'LineWidth',1)
title(leg6,'CMIP6 models','FontName','Arial','FontSize',10)

%%
cLand5_avg11_05kgC = nanmean(cLand5_map11_05kgC,3); 
cLand6_avg11_05kgC = nanmean(cLand6_map11_05kgC,3);

cLand5_avg11_05kgC(cLand5_avg11_05kgC>999) = NaN;
cLand6_avg11_05kgC(cLand6_avg11_05kgC>999) = NaN;

% Observation-derived data
cLand_map_obs = cLand_obs_mean;
% CMIP5 model ensemble mean 
cLand5_avg_M11 = cLand5_avg11_05kgC;
% CMIP6 model ensemble mean 
cLand6_avg_M11 = cLand6_avg11_05kgC;

CMIP5_bias = cLand5_avg_M11 - cLand_map_obs;
CMIP6_bias = cLand6_avg_M11 - cLand_map_obs;

Labels = {'(a) CMIP5','(c) CMIP6'};
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\mycolor_tuaE_bias.mat

%ax = gca
%mycolor_cLand_bias = colormap(ax)
%save ('E:\1_Mycase\3_CMIP56_Cland\4_Version3\1_Figure\2_Figure2_X_GPP_tauE\mycolor_cLand_bias.mat','mycolor_cLand_bias') 
load E:\1_Mycase\3_CMIP56_Cland\4_Version3\1_Figure\2_Figure2_X_GPP_tauE\mycolor_cLand_bias.mat

%figure
%set(gcf,'position',[100 100 490.2,530])
maps2 = tight_subplot(2,1,[-0.08 -0.085],[0.25 -0.018],[0.37 0.09])

cLand5_Bia_M11 = CMIP5_bias;
cLand5_Bia_M11(302:360,:) = [];
cLand5_Bia_M11 = flipud(cLand5_Bia_M11);
raster5_cLand_bia = georasterref('RasterSize',size(cLand5_Bia_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);

cmip5 = maps2(1)
axes(cmip5)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('FLineWidth',1,'FEdgeColor','none')
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3,'DefaultEdgeColor',[0.45 0.45 0.45])
framem('FLineWidth',1)
h = geoshow(cLand5_Bia_M11,raster5_cLand_bia, 'DisplayType','surface','Zdata',zeros(size(cLand5_Bia_M11)),'CData',cLand5_Bia_M11);
colormap(cmip5,mycolor_cLand_bias)
caxis([-20 20])
set(gca,'box','off')
setm(cmip5, 'FFaceColor', [0.68,0.88,0.96])
framem('FLineWidth',1,'FEdgeColor','none')
axis off
text(-2.284,2.0824, '(b) CMIP5','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);

% CMIP6 tuaE bias
cLand6_bia_M11 = CMIP6_bias;
cLand6_bia_M11(302:360,:) = [];
cLand6_bia_M11 = flipud(cLand6_bia_M11);
raster6_cLand_bia = georasterref('RasterSize',size(cLand6_bia_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);   
    
cmip6 = maps2(2)
axes(cmip6)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('FLineWidth',1)
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3,'DefaultEdgeColor',[0.45 0.45 0.45])
framem('FLineWidth',1,'FEdgeColor','none')
h = geoshow(cLand6_bia_M11,raster6_cLand_bia, 'DisplayType','surface','Zdata',zeros(size(cLand6_bia_M11)),'CData',cLand6_bia_M11);
colormap(cmip6,mycolor_cLand_bias)
caxis([-20 20])
set(gca,'box','off')
setm(cmip6, 'FFaceColor', [0.68,0.88,0.96])
framem('FLineWidth',1,'FEdgeColor','none')
axis off
text(-2.284,2.0824, '(e) CMIP6','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);    
h1 = colorbar
h1.Location = 'southoutside'
h1.Position = [0.4753,0.28,0.32714351425943,0.0298];
h1.FontName = 'Arial'
h1.FontSize = 10;
text(0.0183,-2.3264,'Model-data difference in cLand (KgC m^-^2) ',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11)  

     
panel2 = tight_subplot(2,1,[0.15 0.25],[0.32 0.15],[0.475 0.46]) 
axes(panel2(1))
hold on

Obs = 2285;
for i =1:11
    plot([0.5,1],[X5_endT.X(i)-Obs X5_endT.X(i)-Obs] ,'LineWidth',3,'LineStyle','-','color',mycolor5(i,:))
    
end
%H_pa5 = patch([0.75,0.75,1.25,1.25],[avgX5-sdX5,avgX5+sdX5,avgX5+sdX5,avgX5-sdX5],[0.6,0.6,0.6]);
%set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
%avg5 = line([0.65 1.25],[avgX5 avgX5],'color','k','linewidth',2);
set(gca,'XLim',[0.4 1.6],'YLim',[-1500 1500],'LineWidth',1.2,'box','off')
set(gca,'Xcolor','none')
set(gca,'Fontname','Arial','FontSize',10)
set(gca,'YTickLabelMode','auto');
%ylabel('Model bias (PgC)','Fontname','Arial','FontSize',11)
set(gca,'color','none','TickDir','out','TickLength',[0.03 0.035])


% CMIP6
axes(panel2(2)) 
hold on
for i =1:11
    plot([0.5,1],[X6_endT.X(i)-Obs X6_endT.X(i)-Obs] ,'LineWidth',3,'LineStyle','-','color',mycolor6(i,:))
    
end
set(gca,'XLim',[0.4 1.6],'YLim',[-1500 1500],'LineWidth',1.2,'box','off')
set(gca,'Xcolor','none')
set(gca,'Fontname','Arial','FontSize',10)
set(gca,'YTickLabelMode','auto');
%ylabel('Model bias (PgC)','Fontname','Arial','FontSize',11)
set(gca,'color','none','TickDir','out','TickLength',[0.03 0.035])    
    
    

%%
cLand5_map11_05kgC(cLand5_map11_05kgC>10^3) = NaN;
cLand6_map11_05kgC(cLand6_map11_05kgC>10^3) = NaN;
cLand5_zonal_11 = nanmean(cLand5_map11_05kgC,2); cLand5_zonal_11 = squeeze(cLand5_zonal_11); cLand5_zonal_11(302:360,:) = [];
cLand6_zonal_11 = nanmean(cLand6_map11_05kgC,2); cLand6_zonal_11 = squeeze(cLand6_zonal_11); cLand6_zonal_11(302:360,:) = [];

cLand_zonal_obs = nanmean(cLand_obs_mean,2); cLand_zonal_obs(302:360,:) = []; 
cLand_zonal_obs4 = nanmean(cLand_obs_range,2); cLand_zonal_obs4 = squeeze(cLand_zonal_obs4);cLand_zonal_obs4(302:360,:) = [];   
cLand_zonal_OBSmin = min(cLand_zonal_obs4,[],2);
cLand_zonal_OBSmax = max(cLand_zonal_obs4,[],2);
boudary_cLand = [cLand_zonal_OBSmin cLand_zonal_OBSmax] - cLand_zonal_obs;
boudary_cLand(:,1) = -boudary_cLand(:,1); boudary_cLand(isnan(boudary_cLand)) = 0;

%CMIP5 
panel3 = tight_subplot(2,1,[0.041 0.01],[0.32 0.04 ],[0.85 0.008])      
axes(panel3(1))
hold on
lat = 90:-0.5:-60;
cLand5_zonal_11(cLand5_zonal_11<=0) = NaN;
for i=1:11
    CMIP5_lines(i) = plot(cLand5_zonal_11(:,i),lat,'color',[0.67,0.91,1.00],'LineWidth',0.9)
end

cLand5_zonal_avg = nanmean(cLand5_zonal_11,2);
cLand5_zonal_SD = nanstd(cLand5_zonal_11,0,2);

Obs_line = boundedline(cLand_zonal_obs,lat,boudary_cLand,'alpha','transparency',0.3,...
'orientation', 'horiz','cmap',[0.15,0.15,0.15],'nan','remove');
CMIP5_avg = plot(cLand5_zonal_avg,lat,'LineWidth',1.8,'color',[0.02,0.64,0.91])
set(Obs_line,'LineWidth',1.6);
set(gca, 'YLim',[-60 90],'XLim',[0 80]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xticks([0 20 40 60 80]);
xticklabels([]);
plot([0, 80],[0 0],'k--','LineWidth',1) 
text(64,80, '(c)','FontName','Arial','FontSize',11)

plot([43 55],[30 30],'k-','LineWidth',1.8)
text(56,30.2,'Obs','Fontname','Arial','Fontsize',9)
plot([43 55],[15 15],'LineWidth',1.8,'color',[0.02,0.64,0.91])
text(56,15.2,'model','Fontname','Arial','Fontsize',9)


axes(panel3(2))
hold on
lat = 90:-0.5:-60;
cLand6_zonal_11(cLand6_zonal_11<=0) = NaN;
for i=1:11
    CMIP6_lines(i) = plot(cLand6_zonal_11(:,i),lat,'color',[0.99,0.76,0.99],'LineWidth',0.9)
end
cLand6_zonal_avg = nanmean(cLand6_zonal_11,2);
cLand6_zonal_SD = nanstd(cLand6_zonal_11,0,2);
%CMIP6_avg = boundedline(GPP6_zonal_avg,lat,GPP6_zonal_SD,'alpha','transparency',0.2,...
%'orientation', 'horiz','cmap',[1.00,0.07,0.65],'nan','remove');
%set(CMIP6_avg,'LineWidth',1.8);
Obs_line = boundedline(cLand_zonal_obs,lat,boudary_cLand,'alpha','transparency',0.3,...
'orientation', 'horiz','cmap',[0.15,0.15,0.15],'nan','remove');
CMIP6_avg = plot(cLand6_zonal_avg,lat,'LineWidth',1.6,'color',[1.00,0.07,0.65])

set(Obs_line,'LineWidth',1.6);
set(gca, 'YLim',[-60 90],'XLim',[0 80]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
%set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xticks([0 20 40 60 ]);
plot([0, 80],[0 0],'k--','LineWidth',1)
xlabel(['cLand',newline,' (KgC m^-^2)'],'Fontname','Arial','FontSize',11,'Position',[39.3,-86.98,-1])
text(65,80, '(f)','FontName','Arial','FontSize',11)

plot([43 55],[30 30],'k-','LineWidth',1.8)
text(56,30.2,'Obs','Fontname','Arial','Fontsize',9)
plot([43 55],[15 15],'LineWidth',1.8,'color',[1.00,0.07,0.65])
text(56,15.2,'model','Fontname','Arial','Fontsize',9)


length(find(X6_endT.X >= 1516.4 && X6_endT.X <= 2835.2))


X6_endT.X >= 1516.4 & X6_endT.X <= 2835.2 

X5_endT.X >= 1516.4 & X5_endT.X <= 2835.2 

