% Figure2: data-model comparison on GPP
clear;clc;

%%
% load NPP from CMIP5 and CMIP6, 
% time: 2001-2005 ; unit: KgC m-2 yr-1
cd('G:\CMIP5\4_MatData\spatial_data\hist_rcp85')
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
    
clearvars -except sp5_NPP5_11 sp6_NPP5_11
% change modeled NPP to match global map              
NPP5_map11(1:180,1:360,1:11) = NaN;
NPP6_map11(1:180,1:360,1:11) = NaN;
for i =1:11
    i
    
    NPP5_M = sp5_NPP5_11(:,:,i);
    NPP6_M = sp6_NPP5_11(:,:,i);
    
    NPP5_t = NPP5_M';
    NPP6_t = NPP6_M';
    
    NPP5_LR(1:180,1:360) = NaN;
    NPP6_LR(1:180,1:360) = NaN;
    
    NPP5_LR(:,1:180) = NPP5_t(:,181:360);
    NPP5_LR(:,181:360) = NPP5_t(:,1:180);
    NPP6_LR(:,1:180) = NPP6_t(:,181:360);
    NPP6_LR(:,181:360) = NPP6_t(:,1:180);
    
    map_NPP5(1:180,1:360) = NaN;
    map_NPP6(1:180,1:360) = NaN;
    for k=1:180
        map_NPP5(181-k,:) = NPP5_LR(k,:);
        map_NPP6(181-k,:) = NPP6_LR(k,:);
    end
    
    NPP5_map11(:,:,i) = map_NPP5;
    NPP6_map11(:,:,i) = map_NPP6;
      
end
clearvars -except NPP5_map11 NPP6_map11
save E:\1_Mycase\3_CMIP56_Cland\8_CodeFigures_JClimate\MainFigures\M_Figure4\cmip56_NPP.mat

% load observational data
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


%% Comparing at the regional scale
clearvars -except NPP5_map11 NPP6_map11 NPP_obs05
lat05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lat');
lon05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lon');
lat1 = 89.5:-1:-89.5; lat1 = lat1';
lon1 = -179.5:179.5;  lon1 = lon1';

[x05,y05] = meshgrid(lon05,lat05);
[x1,y1] = meshgrid(lon1,lat1);
% regrided CMIP5 and CMIP6 data into 0.5x0.5 resolution
NPP5_sp05(1:360,1:720,1:11) = NaN;
NPP6_sp05(1:360,1:720,1:11) = NaN;
for i=1:11
    NPP5_M = interp2(x1,y1,NPP5_map11(:,:,i),x05,y05,'linear');
    NPP6_M = interp2(x1,y1,NPP6_map11(:,:,i),x05,y05,'linear');
    
    NPP5_sp05(:,:,i) = NPP5_M;
    NPP6_sp05(:,:,i) = NPP6_M;
end

NoDesert_map_obs = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','tau');
NoDesert_map_obs = NoDesert_map_obs(:,:,1)';   % considering 0-1m soil C
NoDesert_map_obs(~isnan(NoDesert_map_obs)) = 1;
mask = NoDesert_map_obs;
% Unit of NPP(both simulated and observed) KgC m-2 yr-1
NPP_obs05KG = NPP_obs05.*10^(-3);
% omit region where NPP < 0.001 KgC m-2 yr-1
NPP5_sp05(NPP5_sp05<=0.001) = nan;
NPP6_sp05(NPP6_sp05<=0.001) = nan;
NPP_obs05KG(NPP_obs05KG<=0.001) = nan;

clearvars -except NPP5_map11 NPP6_map11 NPP_obs05 ...
                  NPP5_sp05 NPP6_sp05 NPP_obs05KG mask
              
NPP_obs05KG = NPP_obs05KG.*mask;
NPP5_sp05 = NPP5_sp05.*mask;
NPP6_sp05 = NPP6_sp05.*mask;            

% calculate across-model mean and standard deviation
NPP_avg11_cmip5 = nanmean(NPP5_sp05,3);
NPP_avg11_cmip6 = nanmean(NPP6_sp05,3);

NPP_sd11_cmip5 = nanstd(NPP5_sp05,0,3);
NPP_sd11_cmip6 = nanstd(NPP6_sp05,0,3);

% observational mean
NPP_obsKgC_avg = nanmean(NPP_obs05KG,3);
% NPP bias
NPP5_bias_glb = NPP_avg11_cmip5 - NPP_obsKgC_avg;
NPP6_bias_glb = NPP_avg11_cmip6 - NPP_obsKgC_avg;

%% Figure
%ax = gca;
%mycolor_GPPsd = colormap(ax)
%save ('E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure2_GPP\mycolor_GPPsd.mat','mycolor_GPPsd') 
load('E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure2_GPP\mycolor_GPPsd.mat');
load E:\1_Mycase\3_CMIP56_Cland\4_Version3\1_Figure\2_Figure2_X_GPP_tauE\mycolor_cLand_bias.mat
figure
%set(gcf,'position',[100 100 499.2,595.2])
set(gcf,'position',[100 100 499.2,630])
%maps = tight_subplot(2,1,[-0.13 -0.085],[0.35 -0.05],[-0.04 0.26])
maps = tight_subplot(2,1,[-0.1 -0.085],[0.35 -0.05],[-0.04 0.26])

% CMIP5 NPP bias
NPP5_bia_M11 = NPP5_bias_glb;
NPP5_bia_M11(302:360,:) = [];
NPP5_bia_M11 = flipud(NPP5_bia_M11);
raster5_NPP_bias = georasterref('RasterSize',size(NPP5_bia_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);

cmip5 = maps(1)
axes(cmip5)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('FLineWidth',1)
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3)
framem('FLineWidth',1)
h = geoshow(NPP5_bia_M11,raster5_NPP_bias, 'DisplayType','surface','Zdata',zeros(size(NPP5_bia_M11)),'CData',NPP5_bia_M11);
colormap(mycolor_cLand_bias)
caxis([-0.6 0.6])
set(gca,'box','off')
setm(cmip5, 'FFaceColor', [1,1,1])
framem('FLineWidth',1,'FEdgeColor','k')
axis off
text(-2.284,2.0824, '(a) CMIP5','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);


% CMIP6 NPP bias
NPP6_bias_M11 = NPP6_bias_glb;
NPP6_bias_M11(302:360,:) = [];
NPP6_bias_M11 = flipud(NPP6_bias_M11);
raster6_NPP_bias = georasterref('RasterSize',size(NPP6_bias_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);   
    
cmip6 = maps(2)
axes(cmip6)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('FLineWidth',1)
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3)
framem('FLineWidth',1)
h = geoshow(NPP6_bias_M11,raster6_NPP_bias, 'DisplayType','surface','Zdata',zeros(size(NPP6_bias_M11)),'CData',NPP6_bias_M11);
colormap(mycolor_cLand_bias)
caxis([-0.6 0.6])
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
text(-0.0033,-2.2663,'model-data difference in NPP (KgC m^-^2 yr^-^1)',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11)    
    
% Zonal mean plots CMIP5
NPP5_sp05(NPP5_sp05>10^3) = NaN;
NPP6_sp05(NPP6_sp05>10^3) = NaN;
NPP5_zonal_11 = nanmean(NPP5_sp05,2); NPP5_zonal_11 = squeeze(NPP5_zonal_11); NPP5_zonal_11(302:360,:) = [];
NPP6_zonal_11 = nanmean(NPP6_sp05,2); NPP6_zonal_11 = squeeze(NPP6_zonal_11); NPP6_zonal_11(302:360,:) = [];

mycolor5 = [255 0 0; 153 51 255; 237 176 33; ...%BCC CanESM2 CCSM4
            0 197 205; 0 205 0; 207 194 124;...   %Had IPSL Miroc
            255 99 71; 65 105 255;...  %mpi NorM
            0 0 0; 158 131 149; 119 136 153]./255;  %BNU, GFDL, mri 

mycolor6 = [255 0 0; 153 51 255; 237 176 33;...    %BCC-CSM2-MR, CanESM5, CESM2
           0 197 205; 0 205 0; 207 194 124;...     %UKESM1-0-LL, IPSL-CM6A-LR,MIROC-ES2L
           255 99 71; 65 105 255;...               %MPI-ESM1-2-LR,NorESM-LM
            0 0 0; 158 131 149; 119 136 153]./255; %ACCESS-ESM1-5, CNRM-ESM2-1  EC-Earth3-Veg  
 
NPP_zonal_obs = nanmean(NPP_obs05KG,2); NPP_zonal_obs = nanmean(NPP_zonal_obs,3); NPP_zonal_obs(302:360,:) = [];  
NPP_zonal5 = nanmean(NPP_obs05KG,2); NPP_zonal5 = squeeze(NPP_zonal5); NPP_zonal5(302:360,:) = []; 
NPP_zonal_min = min(NPP_zonal5,[],2);
NPP_zonal_max = max(NPP_zonal5,[],2);
boudary_NPP = [NPP_zonal_min NPP_zonal_max] - NPP_zonal_obs;
boudary_NPP(:,1) = -boudary_NPP(:,1); boudary_NPP(isnan(boudary_NPP)) = 0;


panel = tight_subplot(2,1,[0.032 0.01],[0.42 0.01],[0.73 0.02])    
axes(panel(1))
hold on
lat = 90:-0.5:-60;

for i=1:11
    CMIP5_lines(i) = plot(NPP5_zonal_11(:,i),lat,'color',[0.67,0.91,1.00],'LineWidth',0.9)
end

NPP5_zonal_avg = nanmean(NPP5_zonal_11,2);
NPP5_zonal_SD = nanstd(NPP5_zonal_11,0,2);
CMIP5_avg = plot(NPP5_zonal_avg,lat,'LineWidth',1.8,'color',[0.02,0.64,0.91])
%CMIP5_avg = boundedline(GPP5_zonal_avg,lat,GPP5_zonal_SD,'alpha','transparency',0.3,...
%'orientation', 'horiz','cmap',[0.30,0.75,0.93],'nan','remove');
%set(CMIP5_avg,'LineWidth',1.8,'color',[0.02,0.64,0.91]);

Obs_line = boundedline(NPP_zonal_obs,lat,boudary_NPP,'alpha','transparency',0.6,...
'orientation', 'horiz','cmap',[0.15,0.15,0.15],'nan','remove');
set(Obs_line,'LineWidth',1.8);
text(1.47,80, '(b)',...
        'FontName','Arial','FontSize',11)
set(gca, 'YLim',[-60 90],'XLim',[0 1.8]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xticks([0 2 4 6])
xticklabels([]);
plot([0, 6],[0 0],'k--','LineWidth',1)

plot([0.8 1.1],[64 64],'k-','LineWidth',1.8)
text(1.2,64.2,'Obs','Fontname','Arial','Fontsize',10)
plot([0.8 1.1],[52 52],'LineWidth',1.8,'color',[0.02,0.64,0.91])
text(1.2,52,'model','Fontname','Arial','Fontsize',10)
    
axes(panel(2))
hold on
lat = 90:-0.5:-60;
for i=1:11
    CMIP6_lines(i) = plot(NPP6_zonal_11(:,i),lat,'color',[0.99,0.76,0.99],'LineWidth',0.9)
end
NPP6_zonal_avg = nanmean(NPP6_zonal_11,2);
NPP6_zonal_SD = nanstd(NPP6_zonal_11,0,2);
CMIP6_avg = plot(NPP6_zonal_avg,lat,'LineWidth',1.8,'color',[1.00,0.07,0.65])
%CMIP6_avg = boundedline(GPP6_zonal_avg,lat,GPP6_zonal_SD,'alpha','transparency',0.2,...
%'orientation', 'horiz','cmap',[1.00,0.07,0.65],'nan','remove');
%set(CMIP6_avg,'LineWidth',1.8);

Obs_line = boundedline(NPP_zonal_obs,lat,boudary_NPP,'alpha','transparency',0.6,...
'orientation', 'horiz','cmap',[0.15,0.15,0.15],'nan','remove');
set(Obs_line,'LineWidth',1.8);
text(1.47,80, '(d)',...
        'FontName','Arial','FontSize',11)
set(gca, 'YLim',[-60 90],'XLim',[0 1.8]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
%set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xticks([0 0.6 1.2 1.8]);
plot([0, 1.8],[0 0],'k--','LineWidth',1)  
xlabel(['NPP',newline,'(KgC m^-^2 yr^-^1)'],'FontName','Arial','FontSize',11)
plot([0.8 1.1],[64 64],'k-','LineWidth',1.8)
text(1.2,64.2,'Obs','Fontname','Arial','Fontsize',10)
plot([0.8 1.1],[52 52],'LineWidth',1.8,'color',[1.00,0.07,0.65])
text(1.2,52,'model','Fontname','Arial','Fontsize',10)

%% circumpolar, non-circumpolar and global terrestrial GPP
clearvars -except NPP5_sp05 NPP6_sp05 NPP_obs05
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\Circumpolar_cSoil\NCSCDv2_Circumpolar_netCDF_05deg\NCSCD05_1mkg.mat')
polar_mask = NCSCDgb_1mKg;
polar_mask(~isnan(polar_mask)) = 1;     % mask for circumpolar region

Noploar_mask = polar_mask;
Noploar_mask(Noploar_mask == 1) = 0;
Noploar_mask(isnan(Noploar_mask)) = 1;
Noploar_mask(Noploar_mask == 0) = NaN;  % mask for non-circumpolar region

% load CMIP5 data
cd('G:\CMIP5\4_MatData\temporal_data\hist_rcp85')
load('NPP_cmip5_tmp.mat')
NaNgf = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
NPPgf_tmp = [NaNgf; NPPgf_tmp];
NaNhad = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
NPPhad_tmp = [NaNhad; NPPhad_tmp];
NPPmri_tmp = [NaN; NPPmri_tmp];
NPP5_tmp(1:155,1:11) = NaN;
NPP5_tmp(:,1) = NPPbcc_tmp(2:156);
NPP5_tmp(:,2) = NPPcan_tmp(2:156);
NPP5_tmp(:,3) = NPPccsm_tmp(2:156);
NPP5_tmp(:,4) = NPPhad_tmp(2:156);
NPP5_tmp(:,5) = NPPipsl_tmp(2:156);
NPP5_tmp(:,6) = NPPmiroc_tmp(2:156);
NPP5_tmp(:,7) = NPPmpi_tmp(2:156);
NPP5_tmp(:,8) = NPPnor_tmp(2:156);
NPP5_tmp(:,9) = NPPbnu_tmp(2:156);
NPP5_tmp(:,10) = NPPgf_tmp(2:156);
NPP5_tmp(:,11) = NPPmri_tmp(2:156);
 
% Load CMIP6 data
cd('H:\CMIP56_Csink\4_MatData\temporal_data\hist_ssp585')
load('2NPP_cmip6_tmp.mat')
NPP6_tmp = [NPPbcc_tmp(2:156),NPPcan_tmp(2:156),NPPcesm_tmp(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   NPPuk_tmp(2:156),NPPipsl_tmp(2:156),NPPmic_tmp(2:156),...          % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   NPPmpi_tmp(2:156),NPPnor_tmp(2:156),...                           % MPI-ESM1-2-LR, NorESM2
   NPPass_tmp(2:156),NPPcnrm_tmp(2:156),NPPec_tmp(2:156)];            % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

Models5 = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'};
Models6 = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'}; 
% global terresrtial GPP (PgC yr-1)    
% models
NPP5_end5_ag = nanmean(NPP5_tmp(151:155,:),1);
NPP6_end5_ag = nanmean(NPP6_tmp(151:155,:),1); 
% observation
area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv');
NPP_obs05_Map_Pg = NPP_obs05.* area05.*10^(-15);    % Unit: PgC yr-1
NPP_obs05_Map_Pg(NPP_obs05_Map_Pg<=0) = NaN;

NPP_obs05_Pg = nansum(NPP_obs05_Map_Pg,1); NPP_obs05_Pg = nansum(NPP_obs05_Pg,2);
NPP_obs05_Pg = squeeze(NPP_obs05_Pg)

% circumpolar 
% models
% CMIP5
NPP5_sp05_map11_Pg = NPP5_sp05.* area05.*10^(-12);  % Unit: PgC yr-1
NPP5_sp05_map11_Pg(NPP5_sp05_map11_Pg<=0) = NaN;
NPP5_sp05_polar_Pg = NPP5_sp05_map11_Pg.* polar_mask;
NPP5_ALLpolar_Pg = nansum(NPP5_sp05_polar_Pg,1); NPP5_ALLpolar_Pg = nansum(NPP5_ALLpolar_Pg,2);
NPP5_ALLpolar_Pg = squeeze(NPP5_ALLpolar_Pg)
% CMIP6
NPP6_sp05_map11_Pg = NPP6_sp05.* area05.*10^(-12);  % Unit: PgC yr-1
NPP6_sp05_map11_Pg(NPP6_sp05_map11_Pg<=0) = NaN;
NPP6_sp05_polar_Pg = NPP6_sp05_map11_Pg.* polar_mask;
NPP6_ALLpolar_Pg = nansum(NPP6_sp05_polar_Pg,1); NPP6_ALLpolar_Pg = nansum(NPP6_ALLpolar_Pg,2);
NPP6_ALLpolar_Pg = squeeze(NPP6_ALLpolar_Pg)
% Observation
NPP_obs05_polar_Pg = NPP_obs05_Map_Pg.* polar_mask;
NPP_ALLpolar_Pg = nansum(NPP_obs05_polar_Pg,1); NPP_ALLpolar_Pg = nansum(NPP_ALLpolar_Pg,2);
NPP_ALLpolar_Pg = squeeze(NPP_ALLpolar_Pg);

% non-circumpolar
% models
NPP5_NonPolarGPP =  NPP5_end5_ag' - NPP5_ALLpolar_Pg;
NPP6_NonPolarGPP =  NPP6_end5_ag' - NPP6_ALLpolar_Pg;
% observation
NPP_obs05_NOpolar_Pg = NPP_obs05_Map_Pg.* Noploar_mask;
NPP_ALL_NOpolar_Pg = nansum(NPP_obs05_NOpolar_Pg,1); NPP_ALL_NOpolar_Pg = nansum(NPP_ALL_NOpolar_Pg,2);
NPP_ALL_NOpolar_Pg = squeeze(NPP_ALL_NOpolar_Pg)

% save GPP5_ALLpolar_Pg GPP6_ALLpolar_Pg GPP5_NonPolarGPP GPP6_NonPolarGPP for the calculation of ecosystem carbon residence time over circumpolar and non-circumpolar regions
save('E:\1_Mycase\3_CMIP56_Cland\8_CodeFigures_JClimate\MainFigures\M_Figure4\NOYpolar_NPP56.mat',...
    'NPP5_ALLpolar_Pg', 'NPP6_ALLpolar_Pg', 'NPP5_NonPolarGPP', 'NPP6_NonPolarGPP','NPP5_end5_ag','NPP6_end5_ag')
    
NL5 = []; NL6 = []; % number of models considered nutrient limitation
NPP5_bar = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','PolarNPP','NonPolarNPP','globalNPP','Nlimitation'});
NPP6_bar = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','PolarNPP','NonPolarNPP','globalNPP','Nlimitation'});

NPP5_bar.Moldes = Models5';
NPP6_bar.Moldes = Models6';
NPP5_bar.PolarNPP = NPP5_ALLpolar_Pg;
NPP6_bar.PolarNPP = NPP6_ALLpolar_Pg;
NPP5_bar.globalNPP = NPP5_end5_ag';
NPP6_bar.globalNPP = NPP6_end5_ag';
NPP5_bar.NonPolarNPP =  NPP5_NonPolarGPP;
NPP6_bar.NonPolarNPP =  NPP6_NonPolarGPP;
NL5 = [0 0 1 0 0 0 0 1 0 0 0]';
NL6 = [0 0 1 1 1 1 1 1  1 0 1]';
NPP5_bar.Nlimitation = NL5;
NPP6_bar.Nlimitation = NL6;

clearvars -except NPP5_bar NPP6_bar NPP_ALLpolar_Pg NPP_ALL_NOpolar_Pg NPP_obs05_Pg
save('E:\1_Mycase\3_CMIP56_Cland\8_CodeFigures_JClimate\MainFigures\M_Figure4\NPP_benchmarkALL.mat',...
    'NPP5_bar','NPP6_bar','NPP_ALLpolar_Pg','NPP_ALL_NOpolar_Pg','NPP_obs05_Pg')

load('E:\1_Mycase\3_CMIP56_Cland\8_CodeFigures_JClimate\MainFigures\M_Figure4\NPP_benchmarkALL.mat')
%%
mycolor5 = [255 0 0; 153 51 255; 237 176 33; ...%BCC CanESM2 CCSM4
            0 197 205; 0 205 0; 207 194 124;...   %Had IPSL Miroc
            255 99 71; 65 105 255;...  %mpi NorM
            0 0 0; 158 131 149; 119 136 153]./255;  %BNU, GFDL, mri 
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


panel2 = tight_subplot(1,3,[0.032 0.08],[0.07 0.70],[0.12 0.04]) 
% Global
avgNPP5 = nanmean(NPP5_bar.globalNPP); avgNPP6 = nanmean(NPP6_bar.globalNPP);
sdNPP5 = nanstd(NPP5_bar.globalNPP); sdNPP6 = nanstd(NPP6_bar.globalNPP);
axes(panel2(1))
hold on
for i = 1:11
    if NPP5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,NPP5_bar.globalNPP(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,NPP5_bar.globalNPP(i),'Marker','>', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if NPP6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,NPP6_bar.globalNPP(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,NPP6_bar.globalNPP(i),'Marker','>', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[avgNPP5-sdNPP5,avgNPP5+sdNPP5,avgNPP5+sdNPP5,avgNPP5-sdNPP5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgNPP6-sdNPP6,avgNPP6+sdNPP6,avgNPP6+sdNPP6,avgNPP6-sdNPP6],[1.00,0.65,0.87]);
line([0.75 1.25],[avgNPP5 avgNPP5],'color','k','linewidth',1.8);
line([1.75 2.25],[avgNPP6 avgNPP6],'color','k','linewidth',1.8);
set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.5)
set(gca,'YLim',[40 100],'XLim',[0.5 2.5],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
x=[0.5 2.5 2.5 0.5];
y=[min(NPP_obs05_Pg) min(NPP_obs05_Pg) max(NPP_obs05_Pg) max(NPP_obs05_Pg)];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
ylabel('NPP (PgC yr^-^1)','Fontname','Arial','FontSize',12)
text(0.55, 94.4,'(e)','Fontname','Arial','FontSize',11)
text(1.2, 36,'Global','Fontname','Arial','FontSize',10)

NPP5_gb_rmse = sqrt(mean((NPP5_bar.globalNPP-mean(NPP_obs05_Pg)).^2))
NPP6_gb_rmse = sqrt(mean((NPP6_bar.globalNPP-mean(NPP_obs05_Pg)).^2))

text(1.37,94.4,'RMSE:','Fontname','Arial','FontSize',8,'color','k')
text(2,94.4,num2str(roundn(NPP5_gb_rmse,-1)),'Fontname','Arial','FontSize',9,'color',[0.00,0.77,0.80])
text(2,89,num2str(roundn(NPP6_gb_rmse,-1)),'Fontname','Arial','FontSize',9,'color',[1.00,0.65,0.87])


% circumpolar
avgNPP5_polar = nanmean(NPP5_bar.PolarNPP); avgNPP6_polar = nanmean(NPP6_bar.PolarNPP);
sdNPP_polar5 = nanstd(NPP5_bar.PolarNPP); sdNPP6_polar = nanstd(NPP6_bar.PolarNPP);
axes(panel2(2))
hold on
for i = 1:11
    if NPP5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,NPP5_bar.PolarNPP(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,NPP5_bar.PolarNPP(i),'Marker','>', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if NPP6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,NPP6_bar.PolarNPP(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,NPP6_bar.PolarNPP(i),'Marker','>', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[avgNPP5_polar-sdNPP_polar5,avgNPP5_polar+sdNPP_polar5,avgNPP5_polar+sdNPP_polar5,avgNPP5_polar-sdNPP_polar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgNPP6_polar-sdNPP6_polar,avgNPP6_polar+sdNPP6_polar,avgNPP6_polar+sdNPP6_polar,avgNPP6_polar-sdNPP6_polar],[1.00,0.65,0.87]);
line([0.75 1.25],[avgNPP5_polar avgNPP5_polar],'color','k','linewidth',1.8);
line([1.75 2.25],[avgNPP6_polar avgNPP6_polar],'color','k','linewidth',1.8);
set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.4)
set(gca,'XLim',[0.5 2.5],'YLim',[2 13],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
x=[0.5 2.5 2.5 0.5];
y=[min(NPP_ALLpolar_Pg) min(NPP_ALLpolar_Pg) max(NPP_ALLpolar_Pg) max(NPP_ALLpolar_Pg)];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
text(0.6, 12,'(f)','Fontname','Arial','FontSize',12)
text(0.9, 1.36,'circumpolar','Fontname','Arial','FontSize',11)

NPP5_polar_rmse = sqrt(mean((NPP5_bar.PolarNPP-mean(NPP_ALLpolar_Pg)).^2))
NPP6_polar_rmse = sqrt(mean((NPP6_bar.PolarNPP-mean(NPP_ALLpolar_Pg)).^2))
text(1.37,12,'RMSE:','Fontname','Arial','FontSize',8,'color','k')
text(2,12,num2str(roundn(NPP5_polar_rmse,-1)),'Fontname','Arial','FontSize',9,'color',[0.00,0.77,0.80])
text(2,11.2,num2str(roundn(NPP6_polar_rmse,-2)),'Fontname','Arial','FontSize',9,'color',[1.00,0.65,0.87])


%non-circumpolar
avgNPP5_NONpolar = nanmean(NPP5_bar.NonPolarNPP); avgGPP6_NONpolar = nanmean(NPP6_bar.NonPolarNPP);
sdNPP_NONpolar5 = nanstd(NPP5_bar.NonPolarNPP); sdGPP6_NONpolar = nanstd(NPP6_bar.NonPolarNPP);
axes(panel2(3))
hold on
for i = 1:11
    if NPP5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,NPP5_bar.NonPolarNPP(i),'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,NPP5_bar.NonPolarNPP(i),'Marker','>', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if NPP6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,NPP6_bar.NonPolarNPP(i),'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,NPP6_bar.NonPolarNPP(i),'Marker','>', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[avgNPP5_NONpolar-sdNPP_NONpolar5,avgNPP5_NONpolar+sdNPP_NONpolar5,avgNPP5_NONpolar+sdNPP_NONpolar5,avgNPP5_NONpolar-sdNPP_NONpolar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgGPP6_NONpolar-sdGPP6_NONpolar,avgGPP6_NONpolar+sdGPP6_NONpolar,avgGPP6_NONpolar+sdGPP6_NONpolar,avgGPP6_NONpolar-sdGPP6_NONpolar],[1.00,0.65,0.87]);
line([0.75 1.25],[avgNPP5_NONpolar avgNPP5_NONpolar],'color','k','linewidth',1.8);
line([1.75 2.25],[avgGPP6_NONpolar avgGPP6_NONpolar],'color','k','linewidth',1.8);
set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.6)
set(gca,'XLim',[0.5 2.5],'YLim',[40 100],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
x=[0.5 2.5 2.5 0.5];
y=[min(NPP_ALL_NOpolar_Pg) min(NPP_ALL_NOpolar_Pg) max(NPP_ALL_NOpolar_Pg) max(NPP_ALL_NOpolar_Pg)];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
text(0.6, 95,'(g)','Fontname','Arial','FontSize',12)
text(0.6, 36.6,'non-circumpolar','Fontname','Arial','FontSize',11)
leg_panel =legend([H_pa5 H_pa6 obs],{'CMIP5','CMIP6','Obs'},'NumColumns',3)
set(leg_panel,'FontName','Arial','FontSize',10,'Position',[0.2848,0.0026,0.5196,0.0318],...
    'color','w','EdgeColor','k')

NPP5_NONpolar_rmse = sqrt(mean((NPP5_bar.NonPolarNPP-mean(NPP_ALL_NOpolar_Pg)).^2))
NPP6_NONpolar_rmse = sqrt(mean((NPP6_bar.NonPolarNPP-mean(NPP_ALL_NOpolar_Pg)).^2))
text(1.37,94.4,'RMSE:','Fontname','Arial','FontSize',8,'color','k')
text(2,94.4,num2str(roundn(NPP5_NONpolar_rmse,-1)),'Fontname','Arial','FontSize',9,'color',[0.00,0.77,0.80])
text(2,89,num2str(roundn(NPP6_NONpolar_rmse,-1)),'Fontname','Arial','FontSize',9,'color',[1.00,0.65,0.87])
