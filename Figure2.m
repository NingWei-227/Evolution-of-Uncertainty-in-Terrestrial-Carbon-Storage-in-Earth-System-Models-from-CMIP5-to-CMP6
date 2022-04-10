% supplementary figure:
% Data-model comparison on cLand over circumpolar and non-cricumpolar regions 
% the circumpolar region was extracted based on NCSCDv2 data
% soil C at 0-1m depth: all models 
% full-depth soil C: CESM2 and NorESM2

% observation-derived data
clear;clc;
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\Biomass\ORNL\Global_Maps_C_Density_2010_1763\cVeg05_KgCm2.mat')
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\HWSD\HWSD_1247\HWSDv2_05.mat')
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\LandGIS\LandGIS05_1m.mat')
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\SoilGrids\SoilGrid05_1mkg.mat')
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\Circumpolar_cSoil\NCSCDv2_Circumpolar_netCDF_05deg\NCSCD05_1mkg.mat')

% soil data
cSoil_data(:,:,1) = HWSD1m_05;
cSoil_data(:,:,2) = LandGIS_1mKg;
cSoil_data(:,:,3) = SoilGrid_1mkg;
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

Noploar_mask = polar_mask;
Noploar_mask(Noploar_mask == 1) = 0;
Noploar_mask(isnan(Noploar_mask)) = 1;
Noploar_mask(Noploar_mask == 0) = NaN;


clearvars -except cSoil_mean cSoil_data ...
                  cVeg_obs05kg cVeg_obs_un05kg ...
                  cLand_obs_mean cLand_obs_range cLand_obs_min  cLand_obs_max...
                  polar_mask Noploar_mask
% omit deserts and other regions where GPP<0.01 kgC m-2 yr-1 based on the criteria of Carvalhais et al., 2013
% tuaE data
%cLand_map_obs = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','tau');
%tuaE_median1m_obs = cLand_map_obs(:,:,1)';   % considering 0-1m soil C
%tuaE_median1m_obs(~isnan(tuaE_median1m_obs)) = 1;
%Mask_tua = tuaE_median1m_obs;
%cSoil_mean = cSoil_mean.* Mask_tua;
%cSoil_data = cSoil_data.* Mask_tua;

%cVeg_obs05kg = cVeg_obs05kg.* Mask_tua;
%cVeg_obs_un05kg = cVeg_obs_un05kg.* Mask_tua;

%cLand_obs_mean = cLand_obs_mean.* Mask_tua;
%cLand_obs_range = cLand_obs_range.* Mask_tua;
%cLand_obs_min  = cLand_obs_min.* Mask_tua;
%cLand_obs_max = cLand_obs_max.* Mask_tua;

%clearvars -except cSoil_mean cSoil_data ...
%                  cVeg_obs05kg cVeg_obs_un05kg ...
%                  cLand_obs_mean cLand_obs_range cLand_obs_min  cLand_obs_max...
%                  polar_mask Mask_tua
              
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
          
clearvars -except cSoil_mean cSoil_data ...
                  cVeg_obs05kg cVeg_obs_un05kg ...
                  cLand_obs_mean cLand_obs_range cLand_obs_min  cLand_obs_max...
                  polar_mask Mask_tua Noploar_mask...
                  cSoil5_BCC cSoil5_CAN cSoil5_CCSM ...
                  cSoil5_HAD cSoil5_IPSL cSoil5_MIROC ...
                  cSoil5_MPI cSoil5_NOR cSoil5_BNU...
                  cSoil5_GF cSoil5_MRI ...
                  cVeg5_BCC cVeg5_CAN cVeg5_CCSM ...
                  cVeg5_HAD cVeg5_IPSL cVeg5_MIROC ...
                  cVeg5_MPI cVeg5_NOR cVeg5_BNU ...
                  cVeg5_GF cVeg5_MRI ...
                  cSoil6_BCC cSoil6_CAN cSoil6_CESM ...
                  cSoil6_UK cSoil6_IPSL cSoil6_MIC ...
                  cSoil6_MPI cSoil6_NOR cSoil6_ASS ...
                  cSoil6_CNRM cSoil6_EC ...
                  cVeg6_BCC cVeg6_CAN cVeg6_CESM ...
                  cVeg6_UK cVeg6_IPSL cVeg6_MIC ...
                  cVeg6_MPI cVeg6_NOR cVeg6_ASS ...
                  cVeg6_CNRM cVeg6_EC 
              
HAD_NaN(1:360,1:180,1:10) = NaN;
cSoil5_HAD = cat(3,HAD_NaN,cSoil5_HAD);
cVeg5_HAD = cat(3,HAD_NaN,cVeg5_HAD);

GF_NaN(1:360,1:180,1:11) = NaN;
cSoil5_GF = cat(3,GF_NaN,cSoil5_GF);
cVeg5_GF = cat(3,GF_NaN,cVeg5_GF);

MRI_NaN(1:360,1:180,1) = NaN;
cSoil5_MRI = cat(3,MRI_NaN,cSoil5_MRI);
cVeg5_MRI = cat(3,MRI_NaN,cVeg5_MRI);              
                                
% cLand (cVeg and cSoil) from 2001 to 2005 
% unit:KgC m-2
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

cLand11_map5_kgC = nanmean(sp5_cLand_kgC,3); cLand11_map5_kgC = squeeze(cLand11_map5_kgC);      % unit: KgC m-2
cLand11_map6_kgC = nanmean(sp6_cLand_kgC,3); cLand11_map6_kgC = squeeze(cLand11_map6_kgC);      % unit: KgC m-2

cLand5_map11_kgC(1:180,1:360,1:11) = NaN;
cLand6_map11_kgC(1:180,1:360,1:11) = NaN;

% convert to match with global maps
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
% convert to 0.5x0.5 degree, only for cLand
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

% load sftlf of 11 models
% CMIP5
cd('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\sftlf\CMIP5')
sftlf5_BCC = ncread('1gre_sftlf_fx_bcc-csm1-1-m_r0i0p0.nc','sftlf');   
sftlf5_CAN = ncread('1gre_sftlf_fx_CanESM2_r0i0p0.nc','sftlf');        
sftlf5_CCSM = ncread('1gre_sftlf_fx_CCSM4_historical_r0i0p0.nc','sftlf');  
sftlf5_HAD = ncread('1gre_sftlf_fx_HadGEM2-ES_r0i0p0.nc','sftlf');            
sftlf5_IPSL = ncread('1gre_sftlf_fx_IPSL-CM5A-MR_historical_r0i0p0.nc','sftlf');    
sftlf5_MIROC = ncread('1gre_sftlf_fx_MIROC-ESM_r0i0p0.nc','sftlf');  
sftlf5_MPI = ncread('1gre_sftlf_fx_MPI-ESM-MR_historical_r0i0p0.nc','sftlf');    
sftlf5_NOR = ncread('1gre_sftlf_fx_NorESM1-M_historical_r0i0p0.nc','sftlf');       
sftlf5_BNU = ncread('1gre_sftlf_fx_BNU-ESM_r0i0p0.nc','sftlf');  
sftlf5_GF = ncread('1gre_sftlf_fx_GFDL-ESM2G_r0i0p0.nc','sftlf');             
sftlf5_MRI = ncread('1gre_sftlf_fx_MRI-ESM1_historical_r0i0p0.nc','sftlf'); 

% CMIP6
% load sftlf of 11 models
cd('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\sftlf\CMIP6')
sftlf6_BCC = ncread('1_1gre_sftlf_fx_BCC-CSM2-MR_hist-resIPO_r1i1p1f1_gn.nc','sftlf');
sftlf6_CAN = ncread('1gre_sftlf_fx_CanESM5_historical_r1i1p1f1_gn.nc','sftlf');
sftlf6_CESM = ncread('1gre_sftlf_fx_CESM2_historical_r1i1p1f1_gn.nc','sftlf');
sftlf6_UK = ncread('1gre_sftlf_fx_UKESM1-0-LL_piControl_r1i1p1f2_gn.nc','sftlf');
sftlf6_IPSL = ncread('1gre_sftlf_fx_IPSL-CM6A-LR_historical_r1i1p1f1_gr.nc','sftlf');
sftlf6_MIC = ncread('1gre_sftlf_fx_MIROC-ES2L_historical_r1i1p1f2_gn.nc','sftlf');
sftlf6_MPI = ncread('1_1gre_sftlf_fx_MPI-ESM1-2-LR_historical_r1i1p1f1_gn.nc','sftlf');
sftlf6_NOR = ncread('1gre_sftlf_fx_NorESM2-LM_historical_r1i1p1f1_gn.nc','sftlf');
sftlf6_ASS = ncread('1_1gre_sftlf_fx_ACCESS-ESM1-5_historical_r1i1p1f1_gn.nc','sftlf');
sftlf6_CNRM = ncread('1gre_sftlf_fx_CNRM-ESM2-1_amip_r1i1p1f2_gr.nc','sftlf');
sftlf6_EC = ncread('1gre_sftlf_fx_EC-Earth3-Veg_historical_r1i1p1f1_gr.nc','sftlf');


% convert unit from KgC m-2 into PgC
sftlf5_NOR = double(sftlf5_NOR);  
area_1gre = ncread('H:\CMIP56\CMIP5_trace\CMIP5_Models\4_cellarea\area.nc','area'); % unit:km2
% CMIP5
cSoil5_BCC_Pg = cSoil5_BCC.* area_1gre.* sftlf5_BCC.*10^6*10^(-12)*10^(-2);
cSoil5_CAN_Pg = cSoil5_CAN.* area_1gre.* sftlf5_CAN.*10^6*10^(-12)*10^(-2);
cSoil5_CCSM_Pg = cSoil5_CCSM.* area_1gre.* sftlf5_CCSM.*10^6*10^(-12)*10^(-2);
cSoil5_HAD_Pg = cSoil5_HAD.* area_1gre.* sftlf5_HAD.*10^6*10^(-12)*10^(-2);
cSoil5_IPSL_Pg = cSoil5_IPSL.* area_1gre.* sftlf5_IPSL.*10^6*10^(-12)*10^(-2);
cSoil5_MIROC_Pg = cSoil5_MIROC.* area_1gre.* sftlf5_MIROC.*10^6*10^(-12)*10^(-2);
cSoil5_MPI_Pg = cSoil5_MPI.* area_1gre.* sftlf5_MPI.*10^6*10^(-12)*10^(-2);
cSoil5_NOR_Pg = cSoil5_NOR.* area_1gre.* sftlf5_NOR.*10^6*10^(-12)*10^(-2);
cSoil5_BNU_Pg = cSoil5_BNU.* area_1gre.* sftlf5_BNU.*10^6*10^(-12)*10^(-2);
cSoil5_GF_Pg = cSoil5_GF.* area_1gre.* sftlf5_GF.*10^6*10^(-12)*10^(-2);
cSoil5_MRI_Pg = cSoil5_MRI.* area_1gre.* sftlf5_MRI.*10^6*10^(-12)*10^(-2);

cVeg5_BCC_Pg = cVeg5_BCC.* area_1gre.* sftlf5_BCC.*10^6*10^(-12)*10^(-2);
cVeg5_CAN_Pg = cVeg5_CAN.* area_1gre.* sftlf5_CAN.*10^6*10^(-12)*10^(-2);
cVeg5_CCSM_Pg = cVeg5_CCSM.* area_1gre.* sftlf5_CCSM.*10^6*10^(-12)*10^(-2);
cVeg5_HAD_Pg = cVeg5_HAD.* area_1gre.* sftlf5_HAD.*10^6*10^(-12)*10^(-2);
cVeg5_IPSL_Pg = cVeg5_IPSL.* area_1gre.* sftlf5_IPSL.*10^6*10^(-12)*10^(-2);
cVeg5_MIROC_Pg = cVeg5_MIROC.* area_1gre.* sftlf5_MIROC.*10^6*10^(-12)*10^(-2);
cVeg5_MPI_Pg = cVeg5_MPI.* area_1gre.* sftlf5_MPI.*10^6*10^(-12)*10^(-2);
cVeg5_NOR_Pg = cVeg5_NOR.* area_1gre.* sftlf5_NOR.*10^6*10^(-12)*10^(-2);
cVeg5_BNU_Pg = cVeg5_BNU.* area_1gre.* sftlf5_BNU.*10^6*10^(-12)*10^(-2);
cVeg5_GF_Pg = cVeg5_GF.* area_1gre.* sftlf5_GF.*10^6*10^(-12)*10^(-2);
cVeg5_MRI_Pg = cVeg5_MRI.* area_1gre.* sftlf5_MRI.*10^6*10^(-12)*10^(-2);


% cVeg and cSoil from 2001 to 2005
% CMIP5
sp5_cVegPg_ed5(:,:,:,1) = cVeg5_BCC_Pg(:,:,152:156);
sp5_cVegPg_ed5(:,:,:,2) = cVeg5_CAN_Pg(:,:,152:156);
sp5_cVegPg_ed5(:,:,:,3) = cVeg5_CCSM_Pg(:,:,152:156);
sp5_cVegPg_ed5(:,:,:,4) = cVeg5_HAD_Pg(:,:,152:156);
sp5_cVegPg_ed5(:,:,:,5) = cVeg5_IPSL_Pg(:,:,152:156);
sp5_cVegPg_ed5(:,:,:,6) = cVeg5_MIROC_Pg(:,:,152:156);
sp5_cVegPg_ed5(:,:,:,7) = cVeg5_MPI_Pg(:,:,152:156);
sp5_cVegPg_ed5(:,:,:,8) = cVeg5_NOR_Pg(:,:,152:156);
sp5_cVegPg_ed5(:,:,:,9) = cVeg5_BNU_Pg(:,:,152:156);
sp5_cVegPg_ed5(:,:,:,10) = cVeg5_GF_Pg(:,:,152:156);
sp5_cVegPg_ed5(:,:,:,11) = cVeg5_MRI_Pg(:,:,152:156);

sp5_cSoilPg_ed5(:,:,:,1) = cSoil5_BCC_Pg(:,:,152:156);
sp5_cSoilPg_ed5(:,:,:,2) = cSoil5_CAN_Pg(:,:,152:156);
sp5_cSoilPg_ed5(:,:,:,3) = cSoil5_CCSM_Pg(:,:,152:156);
sp5_cSoilPg_ed5(:,:,:,4) = cSoil5_HAD_Pg(:,:,152:156);
sp5_cSoilPg_ed5(:,:,:,5) = cSoil5_IPSL_Pg(:,:,152:156);
sp5_cSoilPg_ed5(:,:,:,6) = cSoil5_MIROC_Pg(:,:,152:156);
sp5_cSoilPg_ed5(:,:,:,7) = cSoil5_MPI_Pg(:,:,152:156);
sp5_cSoilPg_ed5(:,:,:,8) = cSoil5_NOR_Pg(:,:,152:156);
sp5_cSoilPg_ed5(:,:,:,9) = cSoil5_BNU_Pg(:,:,152:156);
sp5_cSoilPg_ed5(:,:,:,10) = cSoil5_GF_Pg(:,:,152:156);
sp5_cSoilPg_ed5(:,:,:,11) = cSoil5_MRI_Pg(:,:,152:156);


% CMIP6
cSoil6_BCC_Pg = cSoil6_BCC.*area_1gre.* sftlf6_BCC.*10^6*10^(-12)*10^(-2);
cSoil6_CAN_Pg = cSoil6_CAN.* area_1gre.* sftlf6_CAN.*10^6*10^(-12)*10^(-2);
cSoil6_CESM_Pg = cSoil6_CESM.* area_1gre.* sftlf6_CESM.*10^6*10^(-12)*10^(-2);
cSoil6_UK_Pg = cSoil6_UK.* area_1gre.* sftlf6_UK.*10^6*10^(-12)*10^(-2);
cSoil6_IPSL_Pg = cSoil6_IPSL.* area_1gre.* sftlf6_IPSL.*10^6*10^(-12)*10^(-2);
cSoil6_MIC_Pg = cSoil6_MIC.* area_1gre.* sftlf6_MIC.*10^6*10^(-12)*10^(-2);
cSoil6_MPI_Pg = cSoil6_MPI.* area_1gre.* sftlf6_MPI.*10^6*10^(-12)*10^(-2);
cSoil6_NOR_Pg = cSoil6_NOR.* area_1gre.* sftlf6_NOR.*10^6*10^(-12)*10^(-2);
cSoil6_ASS_Pg = cSoil6_ASS.* area_1gre.* sftlf6_ASS.*10^6*10^(-12)*10^(-2);
cSoil6_CNRM_Pg = cSoil6_CNRM.* area_1gre.* sftlf6_CNRM.*10^6*10^(-12)*10^(-2);
cSoil6_EC_Pg = cSoil6_EC.* area_1gre.* sftlf6_EC.*10^6*10^(-12)*10^(-2);

cVeg6_BCC_Pg = cVeg6_BCC.*area_1gre.* sftlf6_BCC.*10^6*10^(-12)*10^(-2);
cVeg6_CAN_Pg = cVeg6_CAN.* area_1gre.* sftlf6_CAN.*10^6*10^(-12)*10^(-2);
cVeg6_CESM_Pg = cVeg6_CESM.* area_1gre.* sftlf6_CESM.*10^6*10^(-12)*10^(-2);
cVeg6_UK_Pg = cVeg6_UK.* area_1gre.* sftlf6_UK.*10^6*10^(-12)*10^(-2);
cVeg6_IPSL_Pg = cVeg6_IPSL.* area_1gre.* sftlf6_IPSL.*10^6*10^(-12)*10^(-2);
cVeg6_MIC_Pg = cVeg6_MIC.* area_1gre.* sftlf6_MIC.*10^6*10^(-12)*10^(-2);
cVeg6_MPI_Pg = cVeg6_MPI.* area_1gre.* sftlf6_MPI.*10^6*10^(-12)*10^(-2);
cVeg6_NOR_Pg = cVeg6_NOR.* area_1gre.* sftlf6_NOR.*10^6*10^(-12)*10^(-2);
cVeg6_ASS_Pg = cVeg6_ASS.* area_1gre.* sftlf6_ASS.*10^6*10^(-12)*10^(-2);
cVeg6_CNRM_Pg = cVeg6_CNRM.* area_1gre.* sftlf6_CNRM.*10^6*10^(-12)*10^(-2);
cVeg6_EC_Pg = cVeg6_EC.* area_1gre.* sftlf6_EC.*10^6*10^(-12)*10^(-2);


sp6_cVegPg_ed5(:,:,:,1) = cVeg6_BCC_Pg(:,:,152:156);
sp6_cVegPg_ed5(:,:,:,2) = cVeg6_CAN_Pg(:,:,152:156);
sp6_cVegPg_ed5(:,:,:,3) = cVeg6_CESM_Pg(:,:,152:156);
sp6_cVegPg_ed5(:,:,:,4) = cVeg6_UK_Pg(:,:,152:156);
sp6_cVegPg_ed5(:,:,:,5) = cVeg6_IPSL_Pg(:,:,152:156);
sp6_cVegPg_ed5(:,:,:,6) = cVeg6_MIC_Pg(:,:,152:156);
sp6_cVegPg_ed5(:,:,:,7) = cVeg6_MPI_Pg(:,:,152:156);
sp6_cVegPg_ed5(:,:,:,8) = cVeg6_NOR_Pg(:,:,152:156);
sp6_cVegPg_ed5(:,:,:,9) = cVeg6_ASS_Pg(:,:,152:156);
sp6_cVegPg_ed5(:,:,:,10) = cVeg6_CNRM_Pg(:,:,152:156);
sp6_cVegPg_ed5(:,:,:,11) = cVeg6_EC_Pg(:,:,152:156);

sp6_cSoilPg_ed5(:,:,:,1) = cSoil6_BCC_Pg(:,:,152:156);
sp6_cSoilPg_ed5(:,:,:,2) = cSoil6_CAN_Pg(:,:,152:156);
sp6_cSoilPg_ed5(:,:,:,3) = cSoil6_CESM_Pg(:,:,152:156);
sp6_cSoilPg_ed5(:,:,:,4) = cSoil6_UK_Pg(:,:,152:156);
sp6_cSoilPg_ed5(:,:,:,5) = cSoil6_IPSL_Pg(:,:,152:156);
sp6_cSoilPg_ed5(:,:,:,6) = cSoil6_MIC_Pg(:,:,152:156);
sp6_cSoilPg_ed5(:,:,:,7) = cSoil6_MPI_Pg(:,:,152:156);
sp6_cSoilPg_ed5(:,:,:,8) = cSoil6_NOR_Pg(:,:,152:156);
sp6_cSoilPg_ed5(:,:,:,9) = cSoil6_ASS_Pg(:,:,152:156);
sp6_cSoilPg_ed5(:,:,:,10) = cSoil6_CNRM_Pg(:,:,152:156);
sp6_cSoilPg_ed5(:,:,:,11) = cSoil6_EC_Pg(:,:,152:156);

cVeg11_map5_pgC = nanmean(sp5_cVegPg_ed5,3); cVeg11_map5_pgC = squeeze(cVeg11_map5_pgC);        % unit: PgC
cSoil_map5_pgC = nanmean(sp5_cSoilPg_ed5,3); cSoil_map5_pgC = squeeze(cSoil_map5_pgC);          % unit: PgC

cVeg11_map6_pgC = nanmean(sp6_cVegPg_ed5,3); cVeg11_map6_pgC = squeeze(cVeg11_map6_pgC);        % unit: PgC 
cSoil_map6_pgC = nanmean(sp6_cSoilPg_ed5,3); cSoil_map6_pgC = squeeze(cSoil_map6_pgC);          % unit: PgC

 
cVeg5_map11_pgC(1:180,1:360,1:11) = NaN;
cSoil5_map11_pgC(1:180,1:360,1:11) = NaN;

cVeg6_map11_pgC(1:180,1:360,1:11) = NaN;
cSoil6_map11_pgC(1:180,1:360,1:11) = NaN;

% convert to match with circumpolar mask
for i =1:11
    i
    
    Cveg5_avg = cVeg11_map5_pgC(:,:,i);
    Csoil5_avg = cSoil_map5_pgC(:,:,i);
    Cveg6_avg = cVeg11_map6_pgC(:,:,i);
    Csoil6_avg = cSoil_map6_pgC(:,:,i);
    
    Cveg5_t = Cveg5_avg';
    Csoil5_t = Csoil5_avg';
    Cveg6_t = Cveg6_avg';
    Csoil6_t = Csoil6_avg';
    
    Cveg5_LR(1:180,1:360) = NaN;
    Csoil5_LR(1:180,1:360) = NaN;
    Cveg6_LR(1:180,1:360) = NaN;
    Csoil6_LR(1:180,1:360) = NaN;
    
    Cveg5_LR(:,1:180) = Cveg5_t(:,181:360);
    Cveg5_LR(:,181:360) = Cveg5_t(:,1:180);
    Csoil5_LR(:,1:180) = Csoil5_t(:,181:360);
    Csoil5_LR(:,181:360) = Csoil5_t(:,1:180);
    Cveg6_LR(:,1:180) = Cveg6_t(:,181:360);
    Cveg6_LR(:,181:360) = Cveg6_t(:,1:180);
    Csoil6_LR(:,1:180) = Csoil6_t(:,181:360);
    Csoil6_LR(:,181:360) = Csoil6_t(:,1:180);
    
    map_Cveg5(1:180,1:360) = NaN;
    map_Csoil5(1:180,1:360) = NaN;
    map_Cveg6(1:180,1:360) = NaN;
    map_Csoil6(1:180,1:360) = NaN;
    for k=1:180
        map_Cveg5(181-k,:) = Cveg5_LR(k,:);
        map_Csoil5(181-k,:) = Csoil5_LR(k,:);
        map_Cveg6(181-k,:) = Cveg6_LR(k,:);
        map_Csoil6(181-k,:) = Csoil6_LR(k,:);
    end
    
    cVeg5_map11_pgC(:,:,i) = map_Cveg5;
    cSoil5_map11_pgC(:,:,i) = map_Csoil5;
    cVeg6_map11_pgC(:,:,i) = map_Cveg6;
    cSoil6_map11_pgC(:,:,i) = map_Csoil6;
      
end


polar_mask_1gre = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\Circumpolar_cSoil\NCSCDv2_Circumpolar_netCDF_1deg\NCSCDv2_Circumpolar_netCDF_1deg\NCSCDv2_Circumpolar_WGS84_SOCC30_1deg.nc','NCSCDv2');
polar_mask_1gre = double(polar_mask_1gre);
Tglobal(1:360,1:124) = NaN;
polar_mask_1gre = cat(2,polar_mask_1gre,Tglobal); 

polar_maskT = rot90(polar_mask_1gre);
New_mask(1:180,1:360) = NaN;
for k=1:180
    New_mask(181-k,:) = polar_maskT(k,:);
end
New_mask(New_mask<0) = NaN; 
New_mask(~isnan(New_mask)) = 1;
mask_polar_1gre = New_mask; 

New_mask2 = mask_polar_1gre;
New_mask2(New_mask2 == 1) = 0;
New_mask2(isnan(New_mask2)) = 1;
New_mask2(New_mask2 == 0) = NaN;
mask_NOpolar_1gre = New_mask2;

tuaE_median_obs = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','tau');
tuaE_median1m_obs = tuaE_median_obs(:,:,1);

lat05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lat');
lon05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lon');
lat1 = 89.5:-1:-89.5; lat1 = lat1';
lon1 = -179.5:179.5;  lon1 = lon1';
[x05,y05] = meshgrid(lon05,lat05);
[x1,y1] = meshgrid(lon1,lat1);
% regrided original map into 1x1 resolution
tuaE_median1m_1gre = interp2(x05,y05,tuaE_median1m_obs',x1,y1,'linear');

tuaE_median1m_1gre(~isnan(tuaE_median1m_1gre)) = 1;
Mask_tuaE_1gre = tuaE_median1m_1gre;
             
clearvars -except cSoil_mean cSoil_data ...
                  cVeg_obs05kg cVeg_obs_un05kg ...
                  cLand_obs_mean cLand_obs_range cLand_obs_min  cLand_obs_max...
                  polar_mask Mask_tua Noploar_mask...
                  cLand5_map11_05kgC cLand6_map11_05kgC ...
                  cVeg5_map11_pgC cSoil5_map11_pgC cVeg6_map11_pgC cSoil6_map11_pgC ...
                  mask_polar_1gre mask_NOpolar_1gre
              
 %cVeg5_map11_pgC = cVeg5_map11_pgC.* Mask_tuaE_1gre;     % omit deserts and other regions where GPP<0.01 kgC m-2 yr-1 based on the criteria of Carvalhais et al., 2013 
 %cSoil5_map11_pgC = cSoil5_map11_pgC.* Mask_tuaE_1gre;
 %cVeg6_map11_pgC = cVeg6_map11_pgC.* Mask_tuaE_1gre;
 %cSoil6_map11_pgC = cSoil6_map11_pgC.* Mask_tuaE_1gre;
           
% extract circumpolar region and exclude extramely large values
cVeg5_polar_map11 = mask_polar_1gre.*cVeg5_map11_pgC;    cVeg5_polar_map11(cVeg5_polar_map11>10^(20)) = NaN;
cSoil5_polar_map11 = mask_polar_1gre.*cSoil5_map11_pgC;  cSoil5_polar_map11(cSoil5_polar_map11>10^(20)) = NaN;
cVeg6_polar_map11 = mask_polar_1gre.*cVeg6_map11_pgC;    cVeg6_polar_map11(cVeg6_polar_map11>10^(20)) = NaN;
cSoil6_polar_map11 = mask_polar_1gre.*cSoil6_map11_pgC;  cSoil6_polar_map11(cSoil6_polar_map11>10^(20)) = NaN;

cVeg5_polar_avg11 = nansum(cVeg5_polar_map11,1); cVeg5_polar_avg11 = nansum(cVeg5_polar_avg11,2);
cVeg5_polar_avg11 = squeeze(cVeg5_polar_avg11);

cSoil5_polar_avg11 = nansum(cSoil5_polar_map11,1); cSoil5_polar_avg11 = nansum(cSoil5_polar_avg11,2);
cSoil5_polar_avg11 = squeeze(cSoil5_polar_avg11);

cVeg6_polar_avg11 = nansum(cVeg6_polar_map11,1); cVeg6_polar_avg11 = nansum(cVeg6_polar_avg11,2);
cVeg6_polar_avg11 = squeeze(cVeg6_polar_avg11);

cSoil6_polar_avg11 = nansum(cSoil6_polar_map11,1); cSoil6_polar_avg11 = nansum(cSoil6_polar_avg11,2);
cSoil6_polar_avg11 = squeeze(cSoil6_polar_avg11); 

% Benchmarking data from Fan et al., 2020 Earth System Science Data, Table2
tuaE_map_obs = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','tau');
tuaE_median1m_obs = tuaE_map_obs(:,:,1)';   % considering 0-1m soil C
tuaE_median1m_obs(~isnan(tuaE_median1m_obs)) = 1;
Mask_tua = tuaE_median1m_obs;

% Csoil 0-1m
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\HWSD\HWSD_1247\HWSDv2_05.mat')
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\LandGIS\LandGIS05_1m.mat')
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\SoilGrids\SoilGrid05_1mkg.mat')
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\Circumpolar_cSoil\NCSCDv2_Circumpolar_netCDF_05deg\NCSCD05_1mkg.mat')
% unit: KgC m-2  
cSoil_data = cSoil_data.* Mask_tua;
cSoil_data(cSoil_data>320) = NaN;

area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv'); 
Global_soil_obs = cSoil_data.*area05.*10^(-12); Global_soil_obs(Global_soil_obs<= 0) = NaN;
Global_soil_obs = nansum(Global_soil_obs,1); Global_soil_obs = nansum(Global_soil_obs,2); 
Global_soil_obs = squeeze(Global_soil_obs)

Plar_soil_obs = cSoil_data.* polar_mask.*area05.*10^(-12);
NCSCD_pg = NCSCDgb_1mKg.*area05.*10^(-12); 
NCSCD_pg = nansum(NCSCD_pg(:));
Plar_soil_obs = nansum(Plar_soil_obs,1); Plar_soil_obs = nansum(Plar_soil_obs,2);
Plar_soil_obs = squeeze(Plar_soil_obs);
Plar_soil_obs = [Plar_soil_obs; NCSCD_pg]


Global_soil_obs = [2195 2091 1332];
Plar_soil_obs = [796 787 568 567];
NonPlar_soil_obs = [1399 1305 764];

% cVeg
area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv');              
Global_Veg_obs = cVeg_obs05kg.*area05.*10^(-12);
Global_Veg_obs = nansum(Global_Veg_obs(:))
Global_VegUN_obs = cVeg_obs_un05kg.*area05.*10^(-12);
Global_VegUN_obs = nansum(Global_VegUN_obs(:))

Plar_Veg_obs = cVeg_obs05kg.* polar_mask.*area05.*10^(-12);
Plar_Veg_obs = nansum(Plar_Veg_obs(:))
Plar_VegUN_obs = cVeg_obs_un05kg.* polar_mask.*area05.*10^(-12);
Plar_VegUN_obs = nansum(Plar_VegUN_obs(:))

NonPlar_Veg_obs = Global_Veg_obs - Plar_Veg_obs
NonPlar_VegUN_obs = Global_VegUN_obs - Plar_VegUN_obs

% cTotal
cPlar_obs = Plar_soil_obs + Plar_Veg_obs;
cPlar_obs_sd1 = cPlar_obs + Plar_VegUN_obs;
cPlar_obs_sd2 = cPlar_obs - Plar_VegUN_obs;
cPlar_obs = [cPlar_obs_sd1 cPlar_obs_sd2]


cNoPlar_obs = NonPlar_soil_obs + NonPlar_Veg_obs;
cNoPlar_obs_sd1 = cNoPlar_obs + NonPlar_VegUN_obs;
cNoPlar_obs_sd2 = cNoPlar_obs - NonPlar_VegUN_obs;
cNoPlar_obs = [cNoPlar_obs_sd1 cNoPlar_obs_sd2]

cGlobal_obs = Global_soil_obs + Global_Veg_obs;
cGlobal_obs_sd1 = cGlobal_obs + Global_VegUN_obs;
cGlobal_obs_sd2 = cGlobal_obs - Global_VegUN_obs;
cGlobal_obs = [cGlobal_obs_sd1 cGlobal_obs_sd2];



% load global estimations of cVeg and cSoil which were calculated based on native resolution
load E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\4_MatData\cVeg_cSoil_CMIP56_native.mat
cVeg5_org_avg11 = nanmean(cVeg5_org_PG);
cSoil5_org_avg11 = nanmean(cSoil5_org_PG);
cVeg6_org_avg11 = nanmean(cVeg6_org_PG);
cSoil6_org_avg11 = nanmean(cSoil6_org_PG);

cVeg5_NonPlar_avg11 = cVeg5_org_avg11' - cVeg5_polar_avg11;
cSoil5_NonPlar_avg11 = cSoil5_org_avg11' - cSoil5_polar_avg11;

cVeg6_NonPlar_avg11 = cVeg6_org_avg11' - cVeg6_polar_avg11;
cSoil6_NonPlar_avg11 = cSoil6_org_avg11' - cSoil6_polar_avg11;

cPlar5_M11 = cVeg5_polar_avg11 + cSoil5_polar_avg11;
cNonPlar5_M11 = cVeg5_NonPlar_avg11 + cSoil5_NonPlar_avg11;

cPlar6_M11 = cVeg6_polar_avg11 + cSoil6_polar_avg11;
cNonPlar6_M11 = cVeg6_NonPlar_avg11 + cSoil6_NonPlar_avg11;

cLand5_M11 = cVeg5_org_avg11' + cSoil5_org_avg11';
cLand6_M11 = cVeg6_org_avg11' + cSoil6_org_avg11';

clearvars -except cSoil_mean cSoil_data ...
                  cVeg_obs05kg cVeg_obs_un05kg ...
                  cLand_obs_mean cLand_obs_range cLand_obs_min  cLand_obs_max...
                  polar_mask Mask_tua Noploar_mask mask_polar_1gre mask_NOpolar_1gre ...
                  cLand5_map11_05kgC cLand6_map11_05kgC ...
                  cVeg5_polar_avg11 cSoil5_polar_avg11 cPlar5_M11 ...
                  cVeg6_polar_avg11 cSoil6_polar_avg11 cPlar6_M11 ...
                  cVeg5_NonPlar_avg11 cSoil5_NonPlar_avg11 cNonPlar5_M11 ...
                  cVeg6_NonPlar_avg11 cSoil6_NonPlar_avg11 cNonPlar6_M11 ...
                  cLand5_M11 cLand6_M11 ...
                  Plar_Veg_obs Plar_VegUN_obs NonPlar_Veg_obs NonPlar_VegUN_obs cPlar_obs cNoPlar_obs...
                  Global_soil_obs Plar_soil_obs NonPlar_soil_obs
              
%% Figure: carbon storage over the circumpolar regions
% load dataset of colorbar
load H:\CMIP56\Figure_Codes2\Codes_case1\spatial_Codes\MyColor16.mat
load H:\CMIP56\CMIP6_trace\CMIP6_Models\1_Mat_data\5_Models_spatial\MyColor10_X2.mat
Labels = {'(a) Data','(b) CMIP5','(c) CESM6'};

cEco_obs_polar = cLand_obs_mean.* polar_mask;
cEco_map5_polar = cLand5_map11_05kgC.* polar_mask;
cEco_map6_polar = cLand6_map11_05kgC.* polar_mask;
              
sp5_Xpolar_11 = nanmean(cEco_map5_polar,3);
sp6_Xpolar_11 = nanmean(cEco_map6_polar,3);

cmip6FD_Xpolar(:,:,1) = cEco_obs_polar;
cmip6FD_Xpolar(:,:,2) = sp5_Xpolar_11;
cmip6FD_Xpolar(:,:,3) = sp6_Xpolar_11;

figure
set(gcf,'position',[100 100 864 403])
% maps 
maps = tight_subplot(3,1,[-0.02 -0.07],[0.1 0.01],[-0.12 0.4])             
for i=1:3
    
    i
    
    Cpolar_avg = cmip6FD_Xpolar(:,:,i);
    Cpolar_M = Cpolar_avg;
    Cpolar_M(302:360,:) = [];
    Cland_1 = flipud(Cpolar_M);
    raster_C = georasterref('RasterSize',size(Cland_1),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figureID = [1,2,3,4]
    figX56 = maps(figureID(i));
    axes(figX56)
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[45 90])
    %setm(figX56, 'FFaceColor', [0.68,0.88,0.96])
    
    framem('off')
    geoshow('landareas.shp','FaceColor','none','LineWidth',0.3)
    framem('FLineWidth',1,'FEdgeColor',[0,0,0])
    geoshow(Cland_1,raster_C, 'DisplayType','surface','Zdata',zeros(size(Cland_1)),'CData',Cland_1);
    colormap(figX56,mymap10_X2)
    caxis([0 50])
    set(gca,'box','off')
    
    %set(gca,'FaceColor',[0.68,0.88,0.96])
    %set(gca,'FaceColor','blue')
    %set(gca,'LineStyle','none')
    axis off
    colorbar('off')
    hold off
    %text(-0.02352,1.55, Models6{i},'HorizontalAlignment','center',...
    %    'FontName','Arial','FontSize',10);
    text(-3.08,2.07, Labels{figureID(i)},'HorizontalAlignment','left',...
        'FontName','Arial','FontSize',12);
end
h2 = colorbar
h2.Location = 'southoutside'
h2.Position = [0.00927,0.109,0.458,0.0354];
h2.FontName = 'Arial'
h2.FontSize = 10;
text(0.0102,0.1081,'Land C storage (Kg C m^-^2)',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12)
              
% panels
panels = tight_subplot(2,2,[0.1 0.06],[0.04 0.07],[0.54, 0.12])
delete(panels(2));
mycolor5 = [255 0 0; 153 51 255; 237 176 33; ...%BCC CanESM2 CCSM4
            0 197 205; 0 205 0; 207 194 124;...   %Had IPSL Miroc
            255 99 71; 65 105 255;...  %mpi NorM
            0 0 0; 158 131 149; 119 136 153]./255;  %BNU, GFDL, mri 
mycolor6 = [255 0 0; 153 51 255; 237 176 33;...    %BCC-CSM2-MR, CanESM5, CESM2
           0 197 205; 0 205 0; 207 194 124;...     %UKESM1-0-LL, IPSL-CM6A-LR,MIROC-ES2L
           255 99 71; 65 105 255;...               %MPI-ESM1-2-LR,NorESM-LM
            0 0 0; 158 131 149; 119 136 153]./255; %ACCESS-ESM1-5, CNRM-ESM2-1  EC-Earth3-Veg  

leg5_str = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'}              
leg6_str = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'}
    
NL5 = []; NL6 = []; % number of models considered nutrient limitation
X5_polarT = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','Xpolar','VegPolar','SoilPolar','Nlimitation'});
X6_polarT = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','Xpolar','VegPolar','SoilPolar','Nlimitation'});              
X5_polarT.Moldes = leg5_str';
X6_polarT.Moldes = leg6_str';
X5_polarT.Xpolar = cPlar5_M11;
X6_polarT.Xpolar = cPlar6_M11;
X5_polarT.VegPolar = cVeg5_polar_avg11;
X6_polarT.VegPolar = cVeg6_polar_avg11;
X5_polarT.SoilPolar = cSoil5_polar_avg11;
X6_polarT.SoilPolar = cSoil6_polar_avg11;
              
NL5 = [0 0 1 0 0 0 0 1 0 0 0]';
NL6 = [0 0 1 1 1 1 1 1  1 0 1]';

X5_polarT.Nlimitation = NL5;
X6_polarT.Nlimitation = NL6;

% cEcosystem
avgXpolar5 = nanmean(X5_polarT.Xpolar); avgXpolar6 = nanmean(X6_polarT.Xpolar);
sdXpolar5 = nanstd(X5_polarT.Xpolar); sdXpolar6 = nanstd(X6_polarT.Xpolar);

axes(panels(1))
%figure
hold on
for i = 1:11
    if X5_polarT.Nlimitation(i) == 0 
        leg5(i)= plot(1,X5_polarT.Xpolar(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,X5_polarT.Xpolar(i),'Marker','^', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if X6_polarT.Nlimitation(i) == 0 
        leg6(i) = plot(2,X6_polarT.Xpolar(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,X6_polarT.Xpolar(i),'Marker','^', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[avgXpolar5-sdXpolar5,avgXpolar5+sdXpolar5,avgXpolar5+sdXpolar5,avgXpolar5-sdXpolar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgXpolar6-sdXpolar6,avgXpolar6+sdXpolar6,avgXpolar6+sdXpolar6,avgXpolar6-sdXpolar6],[1.00,0.65,0.87]);
line([0.75 1.25],[avgXpolar5 avgXpolar5],'color','k','linewidth',1.8);
line([1.75 2.25],[avgXpolar6 avgXpolar6],'color','k','linewidth',1.8);

set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(gca,'XLim',[0.5 2.5],'YLim',[10 1000],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)

x=[0.5 2.5 2.5 0.5];
y=[min(cPlar_obs) min(cPlar_obs) max(cPlar_obs) max(cPlar_obs)]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.5)
set(gca,'YTickLabelMode','auto');
ylabel('C storage (Pg C)','Fontname','Arial','FontSize',12)
text(0.7, 930,'(d)','Fontname','Arial','FontSize',12)
text(0.715, 1087,'Ecosystem','Fontname','Arial','FontSize',12)

leg_models = legend(leg5,leg5_str,'NumColumns',1)
set(leg_models,'FontName','Arial','FontSize',8,'Position',[0.67,0.519,0.151,0.43],...
    'color','none','EdgeColor','none')

% cPlant
cVegpolar5 = nanmean(X5_polarT.VegPolar); cVegpolar6 = nanmean(X6_polarT.VegPolar);
sd_VegPolar5 = nanstd(X5_polarT.VegPolar); sd_Veg_polar6 = nanstd(X6_polarT.VegPolar);

axes(panels(3))
hold on
for i = 1:11
    if X5_polarT.Nlimitation(i) == 0 
        leg5(i)= plot(1,X5_polarT.VegPolar(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,X5_polarT.VegPolar(i),'Marker','^', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if X6_polarT.Nlimitation(i) == 0 
        leg6(i) = plot(2,X6_polarT.VegPolar(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,X6_polarT.VegPolar(i),'Marker','^', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[cVegpolar5-sd_VegPolar5,cVegpolar5+sd_VegPolar5,cVegpolar5+sd_VegPolar5,cVegpolar5-sd_VegPolar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[cVegpolar6-sd_Veg_polar6,cVegpolar6+sd_Veg_polar6,cVegpolar6+sd_Veg_polar6,cVegpolar6-sd_Veg_polar6],[1.00,0.65,0.87]);
line([0.75 1.25],[cVegpolar5 cVegpolar5],'color','k','linewidth',1.8);
line([1.75 2.25],[cVegpolar6 cVegpolar6],'color','k','linewidth',1.8);

set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(gca,'XLim',[0.5 2.5],'YLim',[0 170],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
%set(gca,'XTick',[0.5:0.5:2.5],'xticklabels',{'', 'CMIP5','','CMIP6',''})
%set(gca,'XTickLabelRotation',15)

x=[0.5 2.5 2.5 0.5];
y=[Plar_Veg_obs-Plar_VegUN_obs Plar_Veg_obs-Plar_VegUN_obs  max(Plar_Veg_obs)+Plar_VegUN_obs max(Plar_Veg_obs)+Plar_VegUN_obs]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.5)
set(gca,'YTickLabelMode','auto');
ylabel('C storage (Pg C)','Fontname','Arial','FontSize',12)
text(0.7, 158,'(e)','Fontname','Arial','FontSize',12)
text(0.715, 184,'Plant','Fontname','Arial','FontSize',12)

leg6_models = legend(leg6,leg6_str,'NumColumns',1)
set(leg6_models,'FontName','Arial','FontSize',8,'Position',[0.81,0.507,0.17,0.441],...
    'color','none','EdgeColor','none')


% cSoil
cSoilpolar5 = nanmean(X5_polarT.SoilPolar); cSoilpolar6 = nanmean(X6_polarT.SoilPolar);
sd_SoilPolar5 = nanstd(X5_polarT.SoilPolar); sd_Soil_polar6 = nanstd(X6_polarT.SoilPolar);

axes(panels(4))
%figure
hold on
for i = 1:11
    if X5_polarT.Nlimitation(i) == 0 
        leg5(i)= plot(1,X5_polarT.SoilPolar(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,X5_polarT.SoilPolar(i),'Marker','^', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if X6_polarT.Nlimitation(i) == 0 
        leg6(i) = plot(2,X6_polarT.SoilPolar(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,X6_polarT.SoilPolar(i),'Marker','^', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[cSoilpolar5-sd_SoilPolar5,cSoilpolar5+sd_SoilPolar5,cSoilpolar5+sd_SoilPolar5,cSoilpolar5-sd_SoilPolar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[cSoilpolar6-sd_Soil_polar6,cSoilpolar6+sd_Soil_polar6,cSoilpolar6+sd_Soil_polar6,cSoilpolar6-sd_Soil_polar6],[1.00,0.65,0.87]);
line([0.75 1.25],[cSoilpolar5 cSoilpolar5],'color','k','linewidth',1.8);
line([1.75 2.25],[cSoilpolar6 cSoilpolar6],'color','k','linewidth',1.8);

set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(gca,'XLim',[0.5 2.5],'YLim',[10 1000],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
%set(gca,'XTick',[0.5:0.5:2.5],'xticklabels',{'', 'CMIP5','','CMIP6',''})
%set(gca,'XTickLabelRotation',15)

x=[0.5 2.5 2.5 0.5];
y=[min(Plar_soil_obs) min(Plar_soil_obs) max(Plar_soil_obs) max(Plar_soil_obs)]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
text(0.7, 930,'(f)','Fontname','Arial','FontSize',12)
text(0.715, 1087,'Soil (0-1m)','Fontname','Arial','FontSize',12)

leg_panel =legend([H_pa5 H_pa6 obs],{'CMIP5','CMIP6','Obs'},'NumColumns',1)
set(leg_panel,'FontName','Arial','FontSize',10,'Position',[0.885,0.263,0.107,0.158],...
    'color','none','EdgeColor','none')

text(-0.165,2319,'CMIP5 models:','Fontname','Arial','FontSize',10)
text(1.884, 2319,'CMIP6 models:','Fontname','Arial','FontSize',10)




%% Figure:  cSoil along full soil depth considered over the circumpolar regions
% CESM2 and NorESM2 simulated soil carbon storage along soil profiles
% we used Full-depth SOM to calculate ecosystem cStorage over circumpolar regions
clearvars -except cSoil_mean cSoil_data ...
                  cVeg_obs05kg cVeg_obs_un05kg ...
                  cLand_obs_mean cLand_obs_range cLand_obs_min  cLand_obs_max...
                  polar_mask Mask_tua Noploar_mask mask_polar_1gre mask_NOpolar_1gre ...
                  cLand5_map11_05kgC cLand6_map11_05kgC ...
                  cVeg5_polar_avg11 cSoil5_polar_avg11 cPlar5_M11 ...
                  cVeg6_polar_avg11 cSoil6_polar_avg11 cPlar6_M11 ...
                  cVeg5_NonPlar_avg11 cSoil5_NonPlar_avg11 cNonPlar5_M11 ...
                  cVeg6_NonPlar_avg11 cSoil6_NonPlar_avg11 cNonPlar6_M11 ...
                  cLand5_M11 cLand6_M11 ...
                  Plar_Veg_obs Plar_VegUN_obs NonPlar_Veg_obs NonPlar_VegUN_obs cPlar_obs cNoPlar_obs...
                  Global_soil_obs Plar_soil_obs NonPlar_soil_obs

cd('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6')
cLit6_CESM = ncread('4_1gre_cLitter_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');
cLit6_NOR = ncread('10_1gre_cLitter_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');
cSOM6_CESM_FD = ncread('4_1gre_cSoil_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cSoil');
cSOM6_NOR_FD = ncread('10_1gre_cSoil_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cSoil');

cd('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6')
cVeg6_CESM = ncread('4_1gre_cVeg_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');
cVeg6_NOR = ncread('10_1gre_cVeg_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');

% 2001-2005 mean cLand
cLand6_CESM_FD = cVeg6_CESM + cLit6_CESM + cSOM6_CESM_FD; cLand6_CESM_FD = nanmean(cLand6_CESM_FD(:,:,152:156),3);
cLand6_NOR_FD = cVeg6_NOR + cLit6_NOR + cSOM6_NOR_FD; cLand6_NOR_FD = nanmean(cLand6_NOR_FD(:,:,152:156),3);

% convert to match with circumpolar mask and extract the circumpolar region
cLand6_CESM_FD = cLand6_CESM_FD';
cLand6_NOR_FD = cLand6_NOR_FD';   
CLand6_CESM_LR(1:180,1:360) = NaN;
CLand6_Nor2_LR(1:180,1:360) = NaN;    
CLand6_CESM_LR(:,1:180) = cLand6_CESM_FD(:,181:360);
CLand6_CESM_LR(:,181:360) = cLand6_CESM_FD(:,1:180);
CLand6_Nor2_LR(:,1:180) = cLand6_NOR_FD(:,181:360);
CLand6_Nor2_LR(:,181:360) = cLand6_NOR_FD(:,1:180);
    
map_cLand_cesmFD(1:180,1:360) = NaN;
map_cLand_norFD(1:180,1:360) = NaN;
for k=1:180     
   map_cLand_cesmFD(181-k,:) = CLand6_CESM_LR(k,:);
   map_cLand_norFD(181-k,:) = CLand6_Nor2_LR(k,:);
end
   
cLand_CESM2_map = map_cLand_cesmFD;    
cLand_Nor2_map = map_cLand_norFD;
cPolar6_CESM_FD = cLand_CESM2_map.* mask_polar_1gre; 
cPolar6_NOR_FD = cLand_Nor2_map.* mask_polar_1gre;

cmip6FD_Xpolar(:,:,1) = cPolar6_CESM_FD;
cmip6FD_Xpolar(:,:,2) = cPolar6_NOR_FD;

load H:\CMIP56\Figure_Codes2\Codes_case1\spatial_Codes\MyColor16.mat
load H:\CMIP56\CMIP6_trace\CMIP6_Models\1_Mat_data\5_Models_spatial\MyColor10_X2.mat
Labels = {'(a) CESM2','(b) NorESM2-LM'};

figure
set(gcf,'position',[100 100 864 266])
% maps 
maps = tight_subplot(2,1,[-0.04 -0.07],[0.15 0.01],[-0.12 0.4])
for i=1:2
    
    i
    
    Cpolar_avg = cmip6FD_Xpolar(:,:,i);
    Cpolar_M = Cpolar_avg;
    Cpolar_M(151:180,:) = [];
    Cland_1 = flipud(Cpolar_M);
    raster_C = georasterref('RasterSize',size(Cland_1),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figureID = [1,2,3,4]
    figX56 = maps(figureID(i));
    axes(figX56)
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[45 90])
    %setm(figX56, 'FFaceColor', [0.68,0.88,0.96])
    
    framem('off')
    geoshow('landareas.shp','FaceColor','none','LineWidth',0.3)
    framem('FLineWidth',1,'FEdgeColor',[0,0,0])
    geoshow(Cland_1,raster_C, 'DisplayType','surface','Zdata',zeros(size(Cland_1)),'CData',Cland_1);
    colormap(figX56,mymap10_X2)
    caxis([0 50])
    set(gca,'box','off')
    
    %set(gca,'FaceColor',[0.68,0.88,0.96])
    %set(gca,'FaceColor','blue')
    %set(gca,'LineStyle','none')
    axis off
    colorbar('off')
    hold off
    %text(-0.02352,1.55, Models6{i},'HorizontalAlignment','center',...
    %    'FontName','Arial','FontSize',10);
    text(-3.08,2.07, Labels{figureID(i)},'HorizontalAlignment','left',...
        'FontName','Arial','FontSize',11);
end
h2 = colorbar
h2.Location = 'southoutside'
h2.Position = [0.0278,0.15663,0.42313,0.0563];
h2.FontName = 'Arial'
h2.FontSize = 11;
text(0.0102,0.1081,'Land C storage (Kg C m^-^2)',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12)


% CESM2 and NorESM2 simulated soil carbon storage along soil profiles
% Benchmarking data from Fan et al., 2020 Earth System Science Data, Table2
% cVeg
area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv');  
Plar_Veg_obs = cVeg_obs05kg.* polar_mask.*area05.*10^(-12);
Plar_Veg_obs = nansum(Plar_Veg_obs(:))
Plar_VegUN_obs = cVeg_obs_un05kg.* polar_mask.*area05.*10^(-12);
Plar_VegUN_obs = nansum(Plar_VegUN_obs(:))

% Csoil full depth
NonPlar_soilFD_obs = [2131 2944 2606];
Plar_soilFD_obs = [1020 1326 1766];
Global_soilFD_obs = [3152 4269 4372];

cPlar_FD_obs = Plar_soilFD_obs + Plar_Veg_obs;
cPlar_FD_obs_SD1 = cPlar_FD_obs + Plar_VegUN_obs;
cPlar_FD_obs_SD2 = cPlar_FD_obs - Plar_VegUN_obs;
cPlar_FD_obs = [cPlar_FD_obs_SD1 cPlar_FD_obs_SD2];


% panels
Panel2 = tight_subplot(1,2,[0.1 0.06],[0.2 0.15],[0.54, 0.12])
% Ecosystem
axes(Panel2(1))
hold on
cPolar_CESM = 1006.3;
cPolar_NOR = 1480.8;
bar(1,cPolar_CESM,'FaceColor',[0.94,0.94,0.94],'EdgeColor',[0.6,0.6,0.6],'LineWidth',1.5);
bar(2,cPolar_NOR,'FaceColor',[0.73,0.91,0.98],'EdgeColor',[0.22,0.75,0.92],'LineWidth',1.5);
bar(3,nanmean(cPlar_FD_obs),'FaceColor',[0.61,0.82,0.73],'EdgeColor',[0.01,0.27,0.10],'LineWidth',1.5);
errorbar(3,nanmean(cPlar_FD_obs),max(cPlar_FD_obs)-nanmean(cPlar_FD_obs),nanmean(cPlar_FD_obs)-min(cPlar_FD_obs),'.','Color','r','LineWidth',1.2 );

set(gca,'linewidth',1,'box','on')
set(gca,'XLim',[0 4],'YLim',[0 2000]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'YTickLabelMode','auto');
%set(gca,'XTickLabelMode','auto');
text(0.2187, 2150,'Ecosystem','Fontname','Arial','FontSize',12)
text(0.2, 1900,'(c)','Fontname','Arial','FontSize',12)
ylabel('Land C storage (Pg C)','Fontname','Arial','FontSize',12)

% Soil
axes(Panel2(2))
hold on
cSoil_Polar_CESM = 943.7;
cSoil_Polar_NOR = 1431.4;

bar_leg(1) = bar(1,cSoil_Polar_CESM,'FaceColor',[0.94,0.94,0.94],'EdgeColor',[0.6,0.6,0.6],'LineWidth',1.5);
bar_leg(2) = bar(2,cSoil_Polar_NOR,'FaceColor',[0.73,0.91,0.98],'EdgeColor',[0.22,0.75,0.92],'LineWidth',1.5);
bar_leg(3) = bar(3,nanmean(Plar_soilFD_obs),'FaceColor',[0.61,0.82,0.73],'EdgeColor',[0.01,0.27,0.10],'LineWidth',1.5);
errorbar(3,nanmean(Plar_soilFD_obs),max(Plar_soilFD_obs)-nanmean(Plar_soilFD_obs),nanmean(Plar_soilFD_obs)-min(Plar_soilFD_obs),'.','Color','r','LineWidth',1.2 );

set(gca,'linewidth',1,'box','on')
set(gca,'XLim',[0 4],'YLim',[0 2000]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'YTickLabelMode','auto');
%set(gca,'XTickLabelMode','auto');
text(0.2187, 2150,'Soil (Full depth)','Fontname','Arial','FontSize',12)
text(0.2, 1900,'(d)','Fontname','Arial','FontSize',12)

leg_Bars =legend(bar_leg,{'CESM2','NorESM2-LM','Observation'},'NumColumns',3)
set(leg_Bars,'FontName','Arial','FontSize',10,'Position',[0.488,0.0838,0.4231,0.0753],...
    'color','none','EdgeColor','k')



%% Figure: carbon storage over non-circumpolar regions (considering soil C at 0-1m depth)
clearvars -except cSoil_mean cSoil_data ...
                  cVeg_obs05kg cVeg_obs_un05kg ...
                  cLand_obs_mean cLand_obs_range cLand_obs_min  cLand_obs_max...
                  polar_mask Mask_tua Noploar_mask mask_polar_1gre mask_NOpolar_1gre ...
                  cLand5_map11_05kgC cLand6_map11_05kgC ...
                  cVeg5_polar_avg11 cSoil5_polar_avg11 cPlar5_M11 ...
                  cVeg6_polar_avg11 cSoil6_polar_avg11 cPlar6_M11 ...
                  cVeg5_NonPlar_avg11 cSoil5_NonPlar_avg11 cNonPlar5_M11 ...
                  cVeg6_NonPlar_avg11 cSoil6_NonPlar_avg11 cNonPlar6_M11 ...
                  cLand5_M11 cLand6_M11 ...
                  Plar_Veg_obs Plar_VegUN_obs NonPlar_Veg_obs NonPlar_VegUN_obs cPlar_obs cNoPlar_obs...
                  Global_soil_obs Plar_soil_obs NonPlar_soil_obs
              
load H:\CMIP56\Figure_Codes2\Codes_case1\spatial_Codes\MyColor16.mat
load H:\CMIP56\CMIP6_trace\CMIP6_Models\1_Mat_data\5_Models_spatial\MyColor10_X2.mat
Labels = {'(a) Data','(b) CMIP5','(c) CMIP6'};

cEco_obs_NOpolar = cLand_obs_mean.* Noploar_mask;
cEco_map5_NOpolar = cLand5_map11_05kgC.* Noploar_mask;
cEco_map6_NOpolar = cLand6_map11_05kgC.* Noploar_mask;

mask_land = cEco_obs_NOpolar;
mask_land(~isnan(mask_land)) = 1;
              
sp5_XNOpolar_11 = nanmean(cEco_map5_NOpolar,3).*mask_land;
sp6_XNOpolar_11 = nanmean(cEco_map6_NOpolar,3).*mask_land;

cmip56_XnonPolar(:,:,1) = cEco_obs_NOpolar;
cmip56_XnonPolar(:,:,2) = sp5_XNOpolar_11;
cmip56_XnonPolar(:,:,3) = sp6_XNOpolar_11;

figure
set(gcf,'position',[100 100 832 688])
% maps 
maps = tight_subplot(3,1,[-0.06 -0.07],[0.06 0],[-0.05 0.4])
for i=1:3
    
    i
    
    CnonPolar_avg = cmip56_XnonPolar(:,:,i);
    CnonPolar_avg(CnonPolar_avg>=9999) = NaN;
    
    CnonPolar_M = CnonPolar_avg;
    CnonPolar_M(302:360,:) = [];
    Cland_1 = flipud(CnonPolar_M);
    raster_C = georasterref('RasterSize',size(Cland_1),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figureID = [1,2,3,4]
    figX56 = maps(figureID(i));
    axes(figX56)
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 70])
    setm(figX56, 'FFaceColor', [0.68,0.88,0.96])
    
    framem('off')
    geoshow('landareas.shp','FaceColor','none','LineWidth',0.3,'DefaultEdgeColor',[0.4 0.4 0.4])
    framem('FLineWidth',1,'FEdgeColor','none')
    geoshow(Cland_1,raster_C, 'DisplayType','surface','Zdata',zeros(size(Cland_1)),'CData',Cland_1);
    colormap(figX56,mymap10_X2)
    caxis([0 50])
    set(gca,'box','off')
    
    %set(gca,'FaceColor',[0.68,0.88,0.96])
    %set(gca,'FaceColor','blue')
    %set(gca,'LineStyle','none')
    axis off
    colorbar('off')
    hold off
    %text(-0.02352,1.55, Models6{i},'HorizontalAlignment','center',...
    %    'FontName','Arial','FontSize',10);
    text(-3.068,-1, Labels{figureID(i)},'HorizontalAlignment','left',...
        'FontName','Arial','FontSize',12);
end
h2 = colorbar
h2.Location = 'southoutside'
h2.Position = [0.02212,0.08605,0.5067,0.0268];
h2.FontName = 'Arial'
h2.FontSize = 12;
text(-0.0253,-2.0932,'Land C storage (Kg C m^-^2)',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',13)

% three panels 
panel_56 = tight_subplot(3,1,[0.07 0.06],[0.15 0.04],[0.61 0.24 ])
mycolor5 = [255 0 0; 153 51 255; 237 176 33; ...%BCC CanESM2 CCSM4
            0 197 205; 0 205 0; 207 194 124;...   %Had IPSL Miroc
            255 99 71; 65 105 255;...  %mpi NorM
            0 0 0; 158 131 149; 119 136 153]./255;  %BNU, GFDL, mri 
mycolor6 = [255 0 0; 153 51 255; 237 176 33;...    %BCC-CSM2-MR, CanESM5, CESM2
           0 197 205; 0 205 0; 207 194 124;...     %UKESM1-0-LL, IPSL-CM6A-LR,MIROC-ES2L
           255 99 71; 65 105 255;...               %MPI-ESM1-2-LR,NorESM-LM
            0 0 0; 158 131 149; 119 136 153]./255; %ACCESS-ESM1-5, CNRM-ESM2-1  EC-Earth3-Veg  

leg5_str = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'}              
leg6_str = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'}
        
NL5 = []; NL6 = []; % number of models considered nutrient limitation
X5_NoPolarT = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','X','Veg','Soil','Nlimitation'});
X6_NoPolarT = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','X','Veg','Soil','Nlimitation'});

X5_NoPolarT.Moldes = leg5_str';
X6_NoPolarT.Moldes = leg6_str';
X5_NoPolarT.X = cNonPlar5_M11;
X6_NoPolarT.X = cNonPlar6_M11;
X5_NoPolarT.Veg = cVeg5_NonPlar_avg11;
X6_NoPolarT.Veg = cVeg6_NonPlar_avg11;
X5_NoPolarT.Soil = cSoil5_NonPlar_avg11;
X6_NoPolarT.Soil = cSoil6_NonPlar_avg11;

NL5 = [0 0 1 0 0 0 0 1 0 0 0]';
NL6 = [0 0 1 1 1 1 1 1  1 0 1]';

X5_NoPolarT.Nlimitation = NL5;
X6_NoPolarT.Nlimitation = NL6;

% cEcosystem
avgX_NonPolar5 = nanmean(X5_NoPolarT.X); avgX_NonPolar6 = nanmean(X6_NoPolarT.X);
sdX_NonPolar5 = nanstd(X5_NoPolarT.X); sdX_NonPolar6 = nanstd(X6_NoPolarT.X);

axes(panel_56(1))
%figure
hold on
for i = 1:11
    if X5_NoPolarT.Nlimitation(i) == 0 
        leg5(i)= plot(1,X5_NoPolarT.X(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,X5_NoPolarT.X(i),'Marker','^', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if X6_NoPolarT.Nlimitation(i) == 0 
        leg6(i) = plot(2,X6_NoPolarT.X(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,X6_NoPolarT.X(i),'Marker','^', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[avgX_NonPolar5-sdX_NonPolar5,avgX_NonPolar5+sdX_NonPolar5,avgX_NonPolar5+sdX_NonPolar5,avgX_NonPolar5-sdX_NonPolar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgX_NonPolar6-sdX_NonPolar6,avgX_NonPolar6+sdX_NonPolar6,avgX_NonPolar6+sdX_NonPolar6,avgX_NonPolar6-sdX_NonPolar6],[1.00,0.65,0.87]);
line([0.75 1.25],[avgX_NonPolar5 avgX_NonPolar5],'color','k','linewidth',1.8);
line([1.75 2.25],[avgX_NonPolar6 avgX_NonPolar6],'color','k','linewidth',1.8);

set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.4)
set(gca,'XLim',[0.5 2.5],'YLim',[500 3400],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',11)

x=[0.5 2.5 2.5 0.5];
y=[min(cNoPlar_obs) min(cNoPlar_obs) max(cNoPlar_obs) max(cNoPlar_obs)]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.3)
set(gca,'YTickLabelMode','auto');
ylabel('C storage (Pg C)','Fontname','Arial','FontSize',12)
%text(0.7, 930,'(e)','Fontname','Arial','FontSize',12)
text(0.715, 3721,'(d) Ecosystem','Fontname','Arial','FontSize',12)

leg_models5 = legend(leg5,leg5_str,'NumColumns',1)
set(leg_models5,'FontName','Arial','FontSize',8,'Position',[0.78,0.646,0.157,0.286],...
    'color','none','EdgeColor','none')
title(leg_models5,'CMIP5 models:','FontSize',10)



% cPlant
cVeg_NoPolar5 = nanmean(X5_NoPolarT.Veg); cVeg_NoPolar6 = nanmean(X6_NoPolarT.Veg);
sd_Veg_NoPolar5 = nanstd(X5_NoPolarT.Veg); sd_Veg_NoPolar6 = nanstd(X6_NoPolarT.Veg);
axes(panel_56(2))
hold on
for i = 1:11
    if X5_NoPolarT.Nlimitation(i) == 0 
        leg5(i)= plot(1,X5_NoPolarT.Veg(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,X5_NoPolarT.Veg(i),'Marker','^', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if X6_NoPolarT.Nlimitation(i) == 0 
        leg6(i) = plot(2,X6_NoPolarT.Veg(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,X6_NoPolarT.Veg(i),'Marker','^', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[cVeg_NoPolar5-sd_Veg_NoPolar5,cVeg_NoPolar5+sd_Veg_NoPolar5,cVeg_NoPolar5+sd_Veg_NoPolar5,cVeg_NoPolar5-sd_Veg_NoPolar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[cVeg_NoPolar6-sd_Veg_NoPolar6,cVeg_NoPolar6+sd_Veg_NoPolar6,cVeg_NoPolar6+sd_Veg_NoPolar6,cVeg_NoPolar6-sd_Veg_NoPolar6],[1.00,0.65,0.87]);
line([0.75 1.25],[cVeg_NoPolar5 cVeg_NoPolar5],'color','k','linewidth',1.8);
line([1.75 2.25],[cVeg_NoPolar6 cVeg_NoPolar6],'color','k','linewidth',1.8);

set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
set(gca,'XLim',[0.5 2.5],'YLim',[0 1000],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
%set(gca,'XTick',[0.5:0.5:2.5],'xticklabels',{'', 'CMIP5','','CMIP6',''})
%set(gca,'XTickLabelRotation',15)

x=[0.5 2.5 2.5 0.5];
y=[NonPlar_Veg_obs-NonPlar_VegUN_obs NonPlar_Veg_obs-NonPlar_VegUN_obs NonPlar_Veg_obs+NonPlar_VegUN_obs NonPlar_Veg_obs+NonPlar_VegUN_obs];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.3)
%text(0.7, 158,'(f)','Fontname','Arial','FontSize',12)
text(0.715, 1080,'(e) Plant','Fontname','Arial','FontSize',12)
ylabel('C storaeg (Pg C)','Fontname','Arial','FontSize',12)


leg_models6 = legend(leg6,leg6_str,'NumColumns',1)
set(leg_models6,'FontName','Arial','FontSize',8,'Position',[0.78,0.274,0.165,0.286],...
    'color','none','EdgeColor','none')
title(leg_models6,'CMIP6 models:','FontSize',10)


% cSoil
cSoil_NoPolar5 = nanmean(X5_NoPolarT.Soil); cSoil_NoPolar6 = nanmean(X6_NoPolarT.Soil);
sd_Soil_NoPolar5 = nanstd(X5_NoPolarT.Soil); sd_Soil_NoPolar6 = nanstd(X6_NoPolarT.Soil);

axes(panel_56(3))
%figure
hold on
for i = 1:11
    if X5_NoPolarT.Nlimitation(i) == 0 
        leg5(i)= plot(1,X5_NoPolarT.Soil(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,X5_NoPolarT.Soil(i),'Marker','^', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if X6_NoPolarT.Nlimitation(i) == 0 
        leg6(i) = plot(2,X6_NoPolarT.Soil(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,X6_NoPolarT.Soil(i),'Marker','^', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[cSoil_NoPolar5-sd_Soil_NoPolar5,cSoil_NoPolar5+sd_Soil_NoPolar5,cSoil_NoPolar5+sd_Soil_NoPolar5,cSoil_NoPolar5-sd_Soil_NoPolar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[cSoil_NoPolar6-sd_Soil_NoPolar6,cSoil_NoPolar6+sd_Soil_NoPolar6,cSoil_NoPolar6+sd_Soil_NoPolar6,cSoil_NoPolar6-sd_Soil_NoPolar6],[1.00,0.65,0.87]);
line([0.75 1.25],[cSoil_NoPolar5 cSoil_NoPolar5],'color','k','linewidth',1.8);
line([1.75 2.25],[cSoil_NoPolar6 cSoil_NoPolar6],'color','k','linewidth',1.8);

set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
set(gca,'XLim',[0.5 2.5],'YLim',[0 3000],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
%set(gca,'XTick',[0.5:0.5:2.5],'xticklabels',{'', 'CMIP5','','CMIP6',''})
%set(gca,'XTickLabelRotation',15)

x=[0.5 2.5 2.5 0.5];
y=[min(NonPlar_soil_obs) min(NonPlar_soil_obs) max(NonPlar_soil_obs) max(NonPlar_soil_obs)]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.3)
%text(0.7, 930,'(g)','Fontname','Arial','FontSize',12)
text(0.715, 3270,'(f) Soil (0-1m)','Fontname','Arial','FontSize',12)
ylabel('C storage (Pg C)','FontName','Arial','FontSize',12)

leg_panel =legend([H_pa5 H_pa6 obs],{'CMIP5','CMIP6','Obs'},'NumColumns',1)
set(leg_panel,'FontName','Arial','FontSize',11,'Position',[0.78,0.157,0.1565,0.081],...
    'color','none','EdgeColor','none')












