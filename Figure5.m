clear;clc;
% benchmark analysis for ecosystem residence time (TuaE)
% Data source of TuaE was calculated based on TuaE = (Ecosystem C storage)/(NPP)
% For consistence, TuaE was calculated in the same way for CMIP56 models in the benchmark analysis

% load land C storage simulated from CMIP5 models
% unit: kgC m-2
cd('G:\CMIP5\4_MatData\spatial_data\hist_rcp85')
load('cLand_cmip5_251.mat')
% calculate 2001-2005 mean
sp5_Cend5(:,:,:,1) = cLand_BCC(:,:,152:156);
sp5_Cend5(:,:,:,2) = cLand_CAN(:,:,152:156);
sp5_Cend5(:,:,:,3) = cLand_CCSM(:,:,152:156);
sp5_Cend5(:,:,:,4) = cLand_HAD(:,:,142:146);

sp5_Cend5(:,:,:,5) = cLand_IPSL(:,:,152:156);
sp5_Cend5(:,:,:,6) = cLand_MIROC(:,:,152:156);
sp5_Cend5(:,:,:,7) = cLand_MPI(:,:,152:156);
sp5_Cend5(:,:,:,8) = cLand_NOR(:,:,152:156);

sp5_Cend5(:,:,:,9) = cLand_BNU(:,:,152:156);
sp5_Cend5(:,:,:,10) = cLand_GF(:,:,141:145);
sp5_Cend5(:,:,:,11) = cLand_MRI(:,:,151:155);

sp5_Cend5_11 = nanmean(sp5_Cend5,3);
sp5_Cend5_11 = squeeze(sp5_Cend5_11);
sp5_Cend5_11(:,1:30,:) = NaN; 

mask = sp5_Cend5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp5_Cend5_11 = sp5_Cend5_11.*mask;

% load NPP data and calculate 2001-2005 mean
% unit: KgC m-2 yr-1
load('NPP_cmip5_251.mat')
sp5_NPP5(:,:,:,1) = nppBCC_3D(:,:,152:156);
sp5_NPP5(:,:,:,2) = nppCAN_3D(:,:,152:156);
sp5_NPP5(:,:,:,3) = nppCCSM_3D(:,:,152:156);
sp5_NPP5(:,:,:,4) = nppHAD_3D(:,:,142:146);

sp5_NPP5(:,:,:,5) = nppIPSL_3D(:,:,152:156);
sp5_NPP5(:,:,:,6) = nppMIROC_3D(:,:,152:156);
sp5_NPP5(:,:,:,7) = nppMPI_3D(:,:,152:156);
sp5_NPP5(:,:,:,8) = nppNOR_3D(:,:,152:156);

sp5_NPP5(:,:,:,9) = nppBNU_3D(:,:,152:156);
sp5_NPP5(:,:,:,10) = nppGF_3D(:,:,141:145);
sp5_NPP5(:,:,:,11) = nppMRI_3D(:,:,151:155);

sp5_NPP5_11 = nanmean(sp5_NPP5,3);
sp5_NPP5_11 = squeeze(sp5_NPP5_11);
sp5_NPP5_11(:,1:30,:) = NaN; 

mask = sp5_NPP5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp5_NPP5_11 = sp5_NPP5_11.*mask;

Models5 = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'};
% omit regions where NPP<= 0.01 KgC m-2 yr-1        
sp5_NPP5_11(sp5_NPP5_11<=0.01)=NaN;
clearvars -except sp5_Cend5_11 sp5_NPP5_11 Models5
% calculate tuaE: tuaE = Cland/NPP, unit: year
numM5 = 11;
for yr = 1:numM5
    
    tuaE_CMIP5(:,:,yr) = sp5_Cend5_11(:,:,yr)./sp5_NPP5_11(:,:,yr);
      
end

clearvars -except sp5_Cend5_11 sp5_NPP5_11 Models5 tuaE_CMIP5
% In CMIP6, CESM2 and NorESM2-LM simulated soil carbon storage along soil depth
% In the benchmark analysis,we only included soil carbon above 1m
% unit: KgC m-2 
cSoil_CESM1m = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\4_1gre_cSoilAbove1m_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cSoilAbove1m');
cSoil_NOR1m = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\10_1gre_cSoilAbove1m_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cSoilAbove1m');

cVeg6_CESM = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\4_1gre_cVeg_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');
cVeg6_NOR = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\10_1gre_cVeg_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');

cLit6_CESM = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\4_1gre_cLitter_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');
cLit6_NOR = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\10_1gre_cLitter_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');

cLand_CESM1m = cVeg6_CESM + cLit6_CESM + cSoil_CESM1m;
cLand_NOR1m = cVeg6_NOR + cLit6_NOR + cSoil_NOR1m;

cd('H:\CMIP56_Csink\4_MatData\spatial_data\hist_ssp585')
load('cLand_cmip6_251.mat')           
sp6_Cend6(:,:,:,1) = cLand_BCC(:,:,152:156);
sp6_Cend6(:,:,:,2) = cLand_CAN(:,:,152:156);
sp6_Cend6(:,:,:,3) = cLand_CESM1m(:,:,152:156);
sp6_Cend6(:,:,:,4) = cLand_UK(:,:,152:156);

sp6_Cend6(:,:,:,5) = cLand_IPSL(:,:,152:156);
sp6_Cend6(:,:,:,6) = cLand_MIC(:,:,152:156);
sp6_Cend6(:,:,:,7) = cLand_MPI(:,:,152:156);
sp6_Cend6(:,:,:,8) = cLand_NOR1m(:,:,152:156);

sp6_Cend6(:,:,:,9) = cLand_ASS(:,:,152:156);
sp6_Cend6(:,:,:,10) = cLand_CNRM(:,:,152:156);
sp6_Cend6(:,:,:,11) = cLand_EC(:,:,152:156);

sp6_Cend5_11 = nanmean(sp6_Cend6,3);
sp6_Cend5_11 = squeeze(sp6_Cend5_11);
sp6_Cend5_11(:,1:30) = NaN;

mask = sp6_Cend5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp6_Cend5_11 = sp6_Cend5_11.*mask;

% load NPP data
% unit: KgC m-2 yr-1
cd('H:\CMIP56_Csink\4_MatData\spatial_data\hist_ssp585')
load('NPP_cmip6_251.mat')
              
sp6_NPP5(:,:,:,1) = nppBCC_3D(:,:,152:156);
sp6_NPP5(:,:,:,2) = nppCAN_3D(:,:,152:156);
sp6_NPP5(:,:,:,3) = nppCESM_3D(:,:,152:156);
sp6_NPP5(:,:,:,4) = nppUK_3D(:,:,152:156);

sp6_NPP5(:,:,:,5) = nppIPSL_3D(:,:,152:156);
sp6_NPP5(:,:,:,6) = nppMIC_3D(:,:,152:156);
sp6_NPP5(:,:,:,7) = nppMPI_3D(:,:,152:156);
sp6_NPP5(:,:,:,8) = nppNOR_3D(:,:,152:156);

sp6_NPP5(:,:,:,9) = nppASS_3D(:,:,152:156);
sp6_NPP5(:,:,:,10) = nppCNRM_3D(:,:,152:156);
sp6_NPP5(:,:,:,11) = nppEC_3D(:,:,152:156);

sp6_NPP5_11 = nanmean(sp6_NPP5,3);
sp6_NPP5_11 = squeeze(sp6_NPP5_11);
sp6_NPP5_11(:,1:30) = NaN;

mask = sp6_NPP5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp6_NPP5_11 = sp6_NPP5_11.*mask;

Models6 = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'};
    
clearvars -except sp5_Cend5_11 sp5_NPP5_11 Models5 tuaE_CMIP5 ...
                  sp6_Cend5_11 sp6_NPP5_11 Models6 
% omit regions where GPP<=0.01 KgC m-2 yr-1             
sp6_NPP5_11(sp6_NPP5_11<=0.01) = NaN;              
num_6M = 11;
tuaE_CMIP6(1:360,1:180,1:11) = NaN;
for yr=1:num_6M 
    tuaE_CMIP6(:,:,yr) = sp6_Cend5_11(:,:,yr)./sp6_NPP5_11(:,:,yr);
    
end

clearvars -except Models5 tuaE_CMIP5 ...
                  Models6 tuaE_CMIP6 

% change tuaE_CMIP5 tuaE_CMIP6 into global map              
tuaE5_map11(1:180,1:360,1:11) = NaN;
tuaE6_map11(1:180,1:360,1:11) = NaN;
for yr =1:11
    yr
    
    tuaE5_M = tuaE_CMIP5(:,:,yr);
    tuaE6_M = tuaE_CMIP6(:,:,yr);
    
    tuaE5_t = tuaE5_M';
    tuaE6_t = tuaE6_M';
    
    tuaE5_LR(1:180,1:360) = NaN;
    tuaE6_LR(1:180,1:360) = NaN;
    
    tuaE5_LR(:,1:180) = tuaE5_t(:,181:360);
    tuaE5_LR(:,181:360) = tuaE5_t(:,1:180);
    tuaE6_LR(:,1:180) = tuaE6_t(:,181:360);
    tuaE6_LR(:,181:360) = tuaE6_t(:,1:180);
    
    map_tuaE5(1:180,1:360) = NaN;
    map_tuaE6(1:180,1:360) = NaN;
    for k=1:180
        map_tuaE5(181-k,:) = tuaE5_LR(k,:);
        map_tuaE6(181-k,:) = tuaE6_LR(k,:);
    end
    
    tuaE5_map11(:,:,yr) = map_tuaE5;
    tuaE6_map11(:,:,yr) = map_tuaE6;
      
end   

% regrid model simulations to match with data 
lat05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lat');
lon05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lon');
lat1 = 89.5:-1:-89.5; lat1 = lat1';
lon1 = -179.5:179.5;  lon1 = lon1';

[x05,y05] = meshgrid(lon05,lat05);
[x1,y1] = meshgrid(lon1,lat1);
% regrided CMIP5 and CMIP6 data into 0.5x0.5 resolution
tuaE5_sp05(1:360,1:720,1:11) = NaN;
tuaE6_sp05(1:360,1:720,1:11) = NaN;
for yr=1:11
    tuaE5_M = interp2(x1,y1,tuaE5_map11(:,:,yr),x05,y05,'linear');
    tuaE6_M = interp2(x1,y1,tuaE6_map11(:,:,yr),x05,y05,'linear');
    
    tuaE5_sp05(:,:,yr) = tuaE5_M;
    tuaE6_sp05(:,:,yr) = tuaE6_M;
end

% observation-based esetimates on ecosystem residence time
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

% Data GIMMS-NPP-GPP
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

% unit: gC m-2 yr-1; 2001-2005 mean 
NPP_obs05(:,:,1) = nanmean(NPP_obs1,3);
NPP_obs05(:,:,2) = nanmean(NPP_obs2,3);
NPP_obs05(:,:,3) = NPP_obs3;
% conver unit into KgC m-2 yr-1
NPP_obs05KG = NPP_obs05.*10^(-3);

% load observation-based estimates on cLand
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
tuaE_mask_fan = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','tau');
tuaE_obs_mean = tuaE_mask_fan(:,:,1)';   % considering 0-1m soil C
tuaE_obs_mean(~isnan(tuaE_obs_mean)) = 1;
Mask_tua = tuaE_obs_mean;

cSoil_mean = cSoil_mean.* Mask_tua;
cSoil_data = cSoil_data.* Mask_tua;

cVeg_obs05kg = cVeg_obs05kg.* Mask_tua;
cVeg_obs_un05kg = cVeg_obs_un05kg.* Mask_tua;

cLand_obs_mean = cLand_obs_mean.* Mask_tua;
cLand_obs_range = cLand_obs_range.* Mask_tua;
cLand_obs_min  = cLand_obs_min.* Mask_tua;
cLand_obs_max = cLand_obs_max.* Mask_tua;

NPP_obs05KG = NPP_obs05KG.* Mask_tua;
tuaE_obs_mean = cLand_obs_mean./nanmean(NPP_obs05KG,3);
tuaE_obs_mean(isinf(tuaE_obs_mean)) = NaN;

tuaE_obs_min(1:360,1:720) = NaN;
tuaE_obs_max(1:360,1:720) = NaN;
for i=1:360
  for j=1:720
      tuaE_obs_min(i,j) = nanmin(cLand_obs_range(i,j,:))./nanmax(NPP_obs05KG(i,j,:));
      tuaE_obs_max(i,j) = nanmax(cLand_obs_range(i,j,:))./nanmin(NPP_obs05KG(i,j,:));
  end
end 
tuaE_obs_min(isinf(tuaE_obs_min)) = NaN;
tuaE_obs_max(isinf(tuaE_obs_max)) = NaN;


tuaE5_sp05 = tuaE5_sp05.* Mask_tua;
tuaE6_sp05 = tuaE6_sp05.* Mask_tua;

clearvars -except Models5 tuaE5_sp05...
                  Models6 tuaE6_sp05...
                  tuaE_obs_min tuaE_obs_max tuaE_obs_mean
save E:\1_Mycase\3_CMIP56_Cland\8_CodeFigures_JClimate\MainFigures\M_Figure5\cmip56_Obs_tuaE.mat           
              
% calculate across-model mean and standard deviation
tuaE5_sp05(tuaE5_sp05==Inf) = NaN; tuaE5_sp05(tuaE5_sp05<= 0) = NaN; tuaE5_sp05(tuaE5_sp05>10^20) = NaN;tuaE5_sp05(tuaE5_sp05==0) = NaN;
tuaE6_sp05(tuaE6_sp05==Inf) = NaN; tuaE6_sp05(tuaE6_sp05<= 0) = NaN; tuaE6_sp05(tuaE6_sp05>10^20) = NaN;tuaE6_sp05(tuaE6_sp05==0) = NaN;
tuaE_avg11_cmip5 = nanmean(tuaE5_sp05,3);
tuaE_avg11_cmip6 = nanmean(tuaE6_sp05,3);
tuaE_sd11_cmip5 = nanstd(tuaE5_sp05,0,3);
tuaE_sd11_cmip6 = nanstd(tuaE6_sp05,0,3); 

tuaE5_bias = tuaE_avg11_cmip5 - tuaE_obs_mean;
tuaE6_bias = tuaE_avg11_cmip6 - tuaE_obs_mean;


%% Figure
%ax = gca;
%mycolor_tuaE_bias = colormap(ax)
%save ('E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\mycolor_tuaE_bias.mat','mycolor_tuaE_bias') 
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\mycolor_tuaE_bias.mat
figure
set(gcf,'position',[100 100 499.2,630])
maps = tight_subplot(2,1,[-0.1 -0.085],[0.35 -0.05],[-0.04 0.26])

% CMIP5 tuaE bias
tuaE5_Bia_M11 = tuaE5_bias;
tuaE5_Bia_M11(302:360,:) = [];
tuaE5_Bia_M11 = flipud(tuaE5_Bia_M11);
raster5_tua_bia = georasterref('RasterSize',size(tuaE5_Bia_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);

cmip5 = maps(1)
axes(cmip5)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('FLineWidth',1)
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3)
framem('FLineWidth',1)
h = geoshow(tuaE5_Bia_M11,raster5_tua_bia, 'DisplayType','surface','Zdata',zeros(size(tuaE5_Bia_M11)),'CData',tuaE5_Bia_M11);
colormap(mycolor_tuaE_bias)
caxis([-100 20])
set(gca,'box','off')
setm(cmip5, 'FFaceColor', [1,1,1])
framem('FLineWidth',1,'FEdgeColor','k')
axis off
text(-2.284,2.0824, '(a) CMIP5','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);

% CMIP6 tuaE bias
tuaE6_bia_M11 = tuaE6_bias;
tuaE6_bia_M11(302:360,:) = [];
tuaE6_bia_M11 = flipud(tuaE6_bia_M11);
raster6_tuaE_bia = georasterref('RasterSize',size(tuaE6_bia_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);   
    
cmip6 = maps(2)
axes(cmip6)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('FLineWidth',1)
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3)
framem('FLineWidth',1)
h = geoshow(tuaE6_bia_M11,raster6_tuaE_bia, 'DisplayType','surface','Zdata',zeros(size(tuaE6_bia_M11)),'CData',tuaE6_bia_M11);
colormap(mycolor_tuaE_bias)
caxis([-100 20])
set(gca,'box','off')
setm(cmip6, 'FFaceColor', [1,1,1])
framem('FLineWidth',1,'FEdgeColor','k')
axis off
text(-2.284,2.0824, '(c) CMIP6','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);    
h1 = colorbar
h1.Location = 'southoutside'
h1.Position = [0.04485,0.3864,0.609,0.0261];
h1.FontName = 'Arial'
h1.FontSize = 10;
text(-0.0033,-2.2663,'model-data difference in \tau_E (yr) ',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11)  

% Zonal mean plots CMIP5
tuaE5_sp05(tuaE5_sp05>10^8) = NaN;
tuaE6_sp05(tuaE6_sp05>10^8) = NaN;
tuaE5_zonal_11 = nanmean(tuaE5_sp05,2); tuaE5_zonal_11 = squeeze(tuaE5_zonal_11); tuaE5_zonal_11(302:360,:) = [];
tuaE6_zonal_11 = nanmean(tuaE6_sp05,2); tuaE6_zonal_11 = squeeze(tuaE6_zonal_11); tuaE6_zonal_11(302:360,:) = [];

mycolor5 = [255 0 0; 153 51 255; 237 176 33; ...%BCC CanESM2 CCSM4
            0 197 205; 0 205 0; 207 194 124;...   %Had IPSL Miroc
            255 99 71; 65 105 255;...  %mpi NorM
            0 0 0; 158 131 149; 119 136 153]./255;  %BNU, GFDL, mri 

mycolor6 = [255 0 0; 153 51 255; 237 176 33;...    %BCC-CSM2-MR, CanESM5, CESM2
           0 197 205; 0 205 0; 207 194 124;...     %UKESM1-0-LL, IPSL-CM6A-LR,MIROC-ES2L
           255 99 71; 65 105 255;...               %MPI-ESM1-2-LR,NorESM-LM
            0 0 0; 158 131 149; 119 136 153]./255; %ACCESS-ESM1-5, CNRM-ESM2-1  EC-Earth3-Veg  
 
tuaE_zonal_obs = nanmean(tuaE_obs_mean,2); tuaE_zonal_obs(302:360,:) = [];  
tuaE_zonal_min = nanmean(tuaE_obs_min,2); tuaE_zonal_min(302:360,:) = [];  
tuaE_zonal_max = nanmean(tuaE_obs_max,2); tuaE_zonal_max(302:360,:) = [];
boudarymin_max = [tuaE_zonal_min tuaE_zonal_max] - tuaE_zonal_obs;
boudarymin_max(:,1) = -boudarymin_max(:,1); boudarymin_max(isnan(boudarymin_max)) = 0;

panel = tight_subplot(2,1,[0.032 0.01],[0.42 0.01],[0.73 0.02])    
axes(panel(1))
hold on
X_info = gca
X_info.XScale = 'log';
X_info.XAxis.Limits = [1, 1000];
X_info.XAxis.TickValues = [10 100 1000];
lat = 90:-0.5:-60;
for yr=1:11
    CMIP5_lines(yr) = plot(tuaE5_zonal_11(:,yr),lat,'color',[0.67,0.91,1.00],'LineWidth',0.9)
end
tuaE5_zonal_avg = nanmean(tuaE5_zonal_11,2);
CMIP5_avg = plot(tuaE5_zonal_avg,lat,'LineWidth',1.8,'color',[0.02,0.64,0.91])
%tuaE5_zonal_SD = nanstd(tuaE5_zonal_11,0,2);
%CMIP5_avg = boundedline(tuaE5_zonal_avg,lat,tuaE5_zonal_SD,'alpha','transparency',0.6,...
%'orientation', 'horiz','cmap',[0.30,0.75,0.93],'nan','remove');
%set(CMIP5_avg,'LineWidth',1.8,'color',[0.02,0.64,0.91]);
Obs_line = boundedline(tuaE_zonal_obs,lat,boudarymin_max,'alpha','transparency',0.4,...
'orientation', 'horiz','cmap',[0.15,0.15,0.15],'nan','remove');
set(Obs_line,'LineWidth',1.8);
text(1.3,82, '(b)',...
        'FontName','Arial','FontSize',11)
set(gca, 'YLim',[-60 90]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xticklabels([]);
plot([1.5, 3.5],[60 60],'k-','LineWidth',1.8)
text(4.2, 60,'Obs','Fontname','Arial','FontSize',9)
plot([1.5, 3.5],[50 50],'LineWidth',1.8,'color',[0.02,0.64,0.91])
text(4.2, 50,'model','Fontname','Arial','FontSize',9)
plot([1, 1000],[0 0],'k--','LineWidth',1)


%leg5 = legend([CMIP5_avg Obs_line],{'Models','Obs'})
%set(leg5,'color','none','EdgeColor','none','Fontname','Arial','Fontsize',10,...
%    'Position',[0.522,0.38,0.1696,0.0392]) 

axes(panel(2))
hold on
X_info = gca
X_info.XAxis.Scale = 'log';
X_info.XAxis.Limits = [1, 1000];
X_info.XAxis.TickValues = [10 100 1000];
lat = 90:-1:-59;
lat = 90:-0.5:-60;
for yr=1:11
    CMIP6_lines(yr) = plot(tuaE6_zonal_11(:,yr),lat,'color',[0.99,0.76,0.99],'LineWidth',0.9)
end
tuaE6_zonal_avg = nanmean(tuaE6_zonal_11,2);
CMIP6_avg = plot(tuaE6_zonal_avg,lat,'LineWidth',1.8,'color',[1.00,0.07,0.65]);
%tuaE6_zonal_SD = nanstd(tuaE6_zonal_11,0,2);

Obs_line = boundedline(tuaE_zonal_obs,lat,boudarymin_max,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.15,0.15,0.15],'nan','remove');
set(Obs_line,'LineWidth',1.8);
text(1.3,82, '(d)',...
        'FontName','Arial','FontSize',11)
set(gca, 'YLim',[-60 90]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
%set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xlabel('\tau_E (yr)', 'FontName','Arial','FontSize',11)  
plot([1.5, 3.5],[60 60],'k-','LineWidth',1.8)
text(4.2, 60,'Obs','Fontname','Arial','FontSize',9)
plot([1.5, 3.5],[50 50],'LineWidth',1.8,'color',[1.00,0.07,0.65])
text(4.2, 50,'model','Fontname','Arial','FontSize',9)
plot([1, 1000],[0 0],'k--','LineWidth',1)

%% comparison on global terrestrial carbon residence time
clear
clc
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\NOY_circumpolar\NOYpolar_cLand.mat
load E:\1_Mycase\3_CMIP56_Cland\8_CodeFigures_JClimate\MainFigures\M_Figure4\NOYpolar_NPP56.mat

% circumpolar
tuaE5_Polar11 = cPlar5_M11./NPP5_ALLpolar_Pg
tuaE6_Polar11 = cPlar6_M11./NPP6_ALLpolar_Pg

% non-circumpolar
tuaE5_NonPolar11 = cNonPlar5_M11./NPP5_NonPolarGPP
tuaE6_NonPolar11 = cNonPlar6_M11./NPP6_NonPolarGPP

% global
tuaE5_end5_ag = cLand5_M11./NPP5_end5_ag'
tuaE6_end5_ag = cLand6_M11./NPP6_end5_ag'

clearvars -except tuaE5_end5_ag tuaE6_end5_ag tuaE5_Polar11 tuaE6_Polar11 tuaE5_NonPolar11 tuaE6_NonPolar11

leg5_str = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'}

leg6_str = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'}

NL5 = []; NL6 = []; % number of models considered nutrient limitation
tuaE5_bar = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','PolartuaE','NoPolartuaE','GlobaltuaE','Nlimitation'});
tuaE6_bar = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','PolartuaE','NoPolartuaE','GlobaltuaE','Nlimitation'});

tuaE5_bar.Moldes = leg5_str';
tuaE6_bar.Moldes = leg6_str';
tuaE5_bar.PolartuaE = tuaE5_Polar11;
tuaE6_bar.PolartuaE = tuaE6_Polar11;
tuaE5_bar.NoPolartuaE = tuaE5_NonPolar11;
tuaE6_bar.NoPolartuaE = tuaE6_NonPolar11;
tuaE5_bar.GlobaltuaE = tuaE5_end5_ag;
tuaE6_bar.GlobaltuaE = tuaE6_end5_ag;
NL5 = [0 0 1 0 0 0 0 1 0 0 0]';
NL6 = [0 0 1 1 1 1 1 1  1 0 1]';
tuaE5_bar.Nlimitation = NL5;
tuaE6_bar.Nlimitation = NL6;

%% Benchmarking data from Fan et al., 2020 Earth System Science Data, Table2
% cVeg (PgC)
NonPlar_Veg_obs = [358 412 399 393];
Plar_Veg_obs = [49 39 38 42];
Global_Veg_obs = [407 451 437 435];

% Csoil 0-1m (PgC)
NonPlar_soil_obs = [1215 1399 1305 764];
Plar_soil_obs = [510 796 787 568 567];
Global_soil_obs = [1725 2195 2091 1332];

% global land carbon storage was estimated based on random combination (PgC)
cLand_obs = [];
cPlar_obs = [];
cNoPlar_obs = [];
for i = 1:4
    cLand  = Global_Veg_obs(i) + Global_soil_obs ;
    C_Plar = Plar_Veg_obs(i) + Plar_soil_obs;
    C_nonPlar = NonPlar_Veg_obs(i) + NonPlar_soil_obs;
    
    cLand_obs = [cLand_obs cLand];
    cPlar_obs = [cPlar_obs C_Plar];
    cNoPlar_obs = [cNoPlar_obs C_nonPlar];
end

% load NPP data
load E:\1_Mycase\3_CMIP56_Cland\8_CodeFigures_JClimate\MainFigures\M_Figure4\NPP_benchmarkALL.mat

% circumpolar
obs_polar_tuaE = [];
for i=1:length(cPlar_obs)
    obs_polar_tuaE(:,i) = cPlar_obs(i)./NPP_ALLpolar_Pg
end
obs_polar_tuaE_mean = nanmean(obs_polar_tuaE(:))
obs_polar_tuaE_SD = nanstd(obs_polar_tuaE(:))

% non-circumpolar
obs_NOpolar_tuaE = [];
for i=1:length(cNoPlar_obs)
    obs_NOpolar_tuaE(:,i) = cNoPlar_obs(i)./NPP_ALL_NOpolar_Pg   
end
obs_NOpolar_tuaE_mean = nanmean(obs_NOpolar_tuaE(:))
obs_NOpolar_tuaE_SD = nanstd(obs_NOpolar_tuaE(:))

% global terresrtial tuaE
obs_glb_tuaE = [];
for i=1:length(cLand_obs)
    obs_glb_tuaE(:,i) = cLand_obs(i)./NPP_obs05_Pg
end
obs_glb_tuaE_mean = nanmean(obs_glb_tuaE(:))
obs_glb_tuaE_SD = nanstd(obs_glb_tuaE(:))

save('E:\1_Mycase\3_CMIP56_Cland\8_CodeFigures_JClimate\MainFigures\M_Figure5\tauE_benchmark.mat',...
    'tuaE5_bar','tuaE6_bar','obs_polar_tuaE_mean','obs_polar_tuaE_SD','obs_NOpolar_tuaE_mean',...
    'obs_NOpolar_tuaE_SD','obs_glb_tuaE_mean','obs_glb_tuaE_SD')

%% Figure
panel2 = tight_subplot(1,3,[0.032 0.08],[0.07 0.70],[0.12 0.04]) 

mycolor5 = [255 0 0; 153 51 255; 237 176 33; ...%BCC CanESM2 CCSM4
            0 197 205; 0 205 0; 207 194 124;...   %Had IPSL Miroc
            255 99 71; 65 105 255;...  %mpi NorM
            0 0 0; 158 131 149; 119 136 153]./255;  %BNU, GFDL, mri 
mycolor6 = [255 0 0; 153 51 255; 237 176 33;...    %BCC-CSM2-MR, CanESM5, CESM2
           0 197 205; 0 205 0; 207 194 124;...     %UKESM1-0-LL, IPSL-CM6A-LR,MIROC-ES2L
           255 99 71; 65 105 255;...               %MPI-ESM1-2-LR,NorESM-LM
            0 0 0; 158 131 149; 119 136 153]./255; %ACCESS-ESM1-5, CNRM-ESM2-1  EC-Earth3-Veg  

% Global
avgtuaE5 = nanmean(tuaE5_bar.GlobaltuaE); avgtuaE6 = nanmean(tuaE6_bar.GlobaltuaE);
sdtuaE5 = nanstd(tuaE5_bar.GlobaltuaE); sdtuaE6 = nanstd(tuaE6_bar.GlobaltuaE);
axes(panel2(1))
hold on
for i = 1:11
    if tuaE5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,tuaE5_bar.GlobaltuaE(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,tuaE5_bar.GlobaltuaE(i),'Marker','>', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if tuaE6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,tuaE6_bar.GlobaltuaE(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,tuaE6_bar.GlobaltuaE(i),'Marker','>', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end

H_pa5 = patch([0.75,0.75,1.25,1.25],[avgtuaE5-sdtuaE5,avgtuaE5+sdtuaE5,avgtuaE5+sdtuaE5,avgtuaE5-sdtuaE5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgtuaE6-sdtuaE6,avgtuaE6+sdtuaE6,avgtuaE6+sdtuaE6,avgtuaE6-sdtuaE6],[1.00,0.65,0.87]);
line([0.75 1.25],[avgtuaE5 avgtuaE5],'color','k','linewidth',1.8);
line([1.75 2.25],[avgtuaE6 avgtuaE6],'color','k','linewidth',1.8);
set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.4)
x=[0.5 2.5 2.5 0.5];
y=[obs_glb_tuaE_mean-obs_glb_tuaE_SD obs_glb_tuaE_mean-obs_glb_tuaE_SD obs_glb_tuaE_mean+obs_glb_tuaE_SD obs_glb_tuaE_mean+obs_glb_tuaE_SD];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YLim',[10 65],'XLim',[0.5 2.5],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
set(gca,'YTickLabelMode','auto');
ylabel('\tau_E (yr)','Fontname','Arial','FontSize',12)
text(0.6, 60,'(e)','Fontname','Arial','FontSize',12)
text(1.2, 6.83,'Global','Fontname','Arial','FontSize',11)

text(1.37,60,'RMSE:','Fontname','Arial','FontSize',8,'color','k')
text(2,60,'17','Fontname','Arial','FontSize',9,'color',[0.00,0.77,0.80])
text(2,55,'15','Fontname','Arial','FontSize',9,'color',[1.00,0.65,0.87])


% circumpolar
avgTuaE5_polar = nanmean(tuaE5_bar.PolartuaE); avgTuaE6_polar = nanmean(tuaE6_bar.PolartuaE);
sdTuaE_polar5 = nanstd(tuaE5_bar.PolartuaE); sdTuaE6_polar = nanstd(tuaE6_bar.PolartuaE);
axes(panel2(2))
hold on
for i = 1:11
    if tuaE5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,tuaE5_bar.PolartuaE(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,tuaE5_bar.PolartuaE(i),'Marker','>', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if tuaE6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,tuaE6_bar.PolartuaE(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,tuaE6_bar.PolartuaE(i),'Marker','>', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end

H_pa5 = patch([0.75,0.75,1.25,1.25],[avgTuaE5_polar-sdTuaE_polar5,avgTuaE5_polar+sdTuaE_polar5,avgTuaE5_polar+sdTuaE_polar5,avgTuaE5_polar-sdTuaE_polar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgTuaE6_polar-sdTuaE6_polar,avgTuaE6_polar+sdTuaE6_polar,avgTuaE6_polar+sdTuaE6_polar,avgTuaE6_polar-sdTuaE6_polar],[1.00,0.65,0.87]);
line([0.75 1.25],[avgTuaE5_polar avgTuaE5_polar],'color','k','linewidth',1.8);
line([1.75 2.25],[avgTuaE6_polar avgTuaE6_polar],'color','k','linewidth',1.8);

set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.4)
set(gca,'XLim',[0.5 2.5],'YLim',[0 250],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)

x=[0.5 2.5 2.5 0.5];
y=[obs_polar_tuaE_mean-obs_polar_tuaE_SD obs_polar_tuaE_mean-obs_polar_tuaE_SD obs_polar_tuaE_mean+obs_polar_tuaE_SD obs_polar_tuaE_mean+obs_polar_tuaE_SD];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
text(0.6, 230,'(f)','Fontname','Arial','FontSize',12)
text(0.9, -13.8,'circumpolar','Fontname','Arial','FontSize',11)

tuaE5_rmse_polar = sqrt(mean(tuaE5_bar.PolartuaE - obs_polar_tuaE_mean).^2)
tuaE6_rmse_polar = sqrt(mean(tuaE6_bar.PolartuaE - obs_polar_tuaE_mean).^2)
text(1.37,230,'RMSE:','Fontname','Arial','FontSize',8,'color','k')
text(2,230,num2str(round(tuaE5_rmse_polar)),'Fontname','Arial','FontSize',9,'color',[0.00,0.77,0.80])
text(2,210,num2str(round(tuaE6_rmse_polar)),'Fontname','Arial','FontSize',9,'color',[1.00,0.65,0.87])



% non-circumpolar
avgTuaE5_NOpolar = nanmean(tuaE5_bar.NoPolartuaE); avgTuaE6_NOpolar = nanmean(tuaE6_bar.NoPolartuaE);
sdTuaE_NOpolar5 = nanstd(tuaE5_bar.NoPolartuaE); sdTuaE6_NOpolar = nanstd(tuaE6_bar.NoPolartuaE);
axes(panel2(3))
hold on
for i = 1:11
    if tuaE5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,tuaE5_bar.NoPolartuaE(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,tuaE5_bar.NoPolartuaE(i),'Marker','>', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if tuaE6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,tuaE6_bar.NoPolartuaE(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,tuaE6_bar.NoPolartuaE(i),'Marker','>', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end

H_pa5 = patch([0.75,0.75,1.25,1.25],[avgTuaE5_NOpolar-sdTuaE_NOpolar5,avgTuaE5_NOpolar+sdTuaE_NOpolar5,avgTuaE5_NOpolar+sdTuaE_NOpolar5,avgTuaE5_NOpolar-sdTuaE_NOpolar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgTuaE6_NOpolar-sdTuaE6_NOpolar,avgTuaE6_NOpolar+sdTuaE6_NOpolar,avgTuaE6_NOpolar+sdTuaE6_NOpolar,avgTuaE6_NOpolar-sdTuaE6_NOpolar],[1.00,0.65,0.87]);
line([0.75 1.25],[avgTuaE5_NOpolar avgTuaE5_NOpolar],'color','k','linewidth',1.8);
line([1.75 2.25],[avgTuaE6_NOpolar avgTuaE6_NOpolar],'color','k','linewidth',1.8);

set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.4)
set(gca,'XLim',[0.5 2.5],'YLim',[10 60],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)

x=[0.5 2.5 2.5 0.5];
y=[obs_NOpolar_tuaE_mean-obs_NOpolar_tuaE_SD obs_NOpolar_tuaE_mean-obs_NOpolar_tuaE_SD obs_NOpolar_tuaE_mean+obs_NOpolar_tuaE_SD obs_NOpolar_tuaE_mean+obs_NOpolar_tuaE_SD];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
text(0.6, 55,'(g)','Fontname','Arial','FontSize',12)
text(0.6, 7.5,'non-circumpolar','Fontname','Arial','FontSize',11)
leg_panel =legend([H_pa5 H_pa6 obs],{'CMIP5','CMIP6','Obs'},'NumColumns',3)
set(leg_panel,'FontName','Arial','FontSize',10,'Position',[0.2848,0.0026,0.5196,0.0318],...
    'color','w','EdgeColor','k')

tuaE5_rmse_NOpolar = sqrt(mean(tuaE5_bar.NoPolartuaE - obs_NOpolar_tuaE_mean).^2)
tuaE6_rmse_NOpolar = sqrt(mean(tuaE6_bar.NoPolartuaE - obs_NOpolar_tuaE_mean).^2)
text(1.37,55,'RMSE:','Fontname','Arial','FontSize',8,'color','k')
text(2,55,num2str(roundn(tuaE5_rmse_NOpolar,-1)),'Fontname','Arial','FontSize',9,'color',[0.00,0.77,0.80])
text(2,50,num2str(roundn(tuaE6_rmse_NOpolar,-1)),'Fontname','Arial','FontSize',9,'color',[1.00,0.65,0.87])
