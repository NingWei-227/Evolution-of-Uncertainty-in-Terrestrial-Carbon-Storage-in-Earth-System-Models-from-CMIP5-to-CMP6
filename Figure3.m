clear
clc

%% global land C storage simulated from Data. CMIP5 and CMIP6
load E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\4_MatData\cVeg_cSoil_CMIP56_native.mat
cVeg5_org_avg11 = nanmean(cVeg5_org_PG);
cSoil5_org_avg11 = nanmean(cSoil5_org_PG);
cVeg6_org_avg11 = nanmean(cVeg6_org_PG);
cSoil6_org_avg11 = nanmean(cSoil6_org_PG);
% CMIP56:
cLand5_Pg_11M = cVeg5_org_avg11 + cSoil5_org_avg11;
cLand6_Pg_11M = cVeg6_org_avg11 + cSoil6_org_avg11; 
% Data
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\Biomass\ORNL\Global_Maps_C_Density_2010_1763\cVeg05_KgCm2.mat')
area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv');              
Global_Veg_obs = cVeg_obs05kg.*area05.*10^(-12);
Global_Veg_obs = nansum(Global_Veg_obs(:))
Global_VegUN_obs = cVeg_obs_un05kg.*area05.*10^(-12);
Global_VegUN_obs = nansum(Global_VegUN_obs(:))
% Benchmarking data from Fan et al., 2020 Earth System Science Data, Table2
% Csoil 0-1m
Global_soil_obs = [2195 2091 1332];

cGlobal_obs = Global_soil_obs + Global_Veg_obs;
cGlobal_obs_sd1 = cGlobal_obs + Global_VegUN_obs;
cGlobal_obs_sd2 = cGlobal_obs - Global_VegUN_obs;
cGlobal_obs = [cGlobal_obs_sd1 cGlobal_obs cGlobal_obs_sd2];

cLand_obs05_Pg = nanmean(cGlobal_obs)

cLand_bias5 = (cLand5_Pg_11M - cLand_obs05_Pg)./cLand_obs05_Pg
cLand_bias6 = (cLand6_Pg_11M - cLand_obs05_Pg)./cLand_obs05_Pg
cLand_rmse5 = sqrt(mean((cLand5_Pg_11M - cLand_obs05_Pg).^2))
cLand_rmse6 = sqrt(mean((cLand6_Pg_11M - cLand_obs05_Pg).^2))

nanmean(cLand_bias5)
nanmean(cLand_bias6)

clearvars -except cLand_bias5 cLand_bias6 cLand_rmse5 cLand_rmse6 cLand_obs05_Pg

%% Global NPP from data-based estimates
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
area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv'); 
NPP_modis_Pg = NPP_obs1.* area05.*10^(-15);    % Unit: PgC yr-1
NPP_modisPG_gb = nansum(NPP_modis_Pg,1); NPP_modisPG_gb = nansum(NPP_modisPG_gb,2); 
NPP_modisPG_gb = squeeze(NPP_modisPG_gb);
NPP_modisPG_Mgb = nanmean(NPP_modisPG_gb); % 5-yr mean

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
NPP_CRUpg = NPP_obs2.* area05.*10^(-15);    % Unit: PgC yr-1
NPP_CRUpg_gb = nansum(NPP_CRUpg,1);NPP_CRUpg_gb = nansum(NPP_CRUpg_gb,2);
NPP_CRUpg_gb = squeeze(NPP_CRUpg_gb);
NPP_CRUpg_Mgb = nanmean(NPP_CRUpg_gb);  % 5-yr mean

% Data CARDAMOM estimates
% reference: Bloom et al., 2015 PNAS
% Unit: gC m-2 day-1
% resolution: 1 x 1
% decadal mean(2001-2010)
npp_CMOMnc = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\DS_10283_875\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_FL_NPP.nc','Mean');
gpp_CMOMnc = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\DS_10283_875\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_OUTPUTv3\CARDAMOM_2001_2010_FL_GPP.nc','Mean');

npp_CMOM_re = npp_CMOMnc';
gpp_CMOM_re = gpp_CMOMnc';
npp_CMOM1gre(1:180,1:360) = NaN;
gpp_CMOM1gre(1:180,1:360) = NaN;
for i=1:180
    npp_CMOM1gre(i,:) =  npp_CMOM_re(181-i,:);
    gpp_CMOM1gre(i,:) =  gpp_CMOM_re(181-i,:);
end 
% regrided CARDAMOM data into 0.5x0.5 resolution
lat05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lat');
lon05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lon');
lat1 = 89.5:-1:-89.5; lat1 = lat1';
lon1 = -179.5:179.5;  lon1 = lon1';
[x05,y05] = meshgrid(lon05,lat05);
[x1,y1] = meshgrid(lon1,lat1);
npp_CMOM05 = interp2(x1,y1,npp_CMOM1gre,x05,y05,'linear');
gpp_CMOM05 = interp2(x1,y1,gpp_CMOM1gre,x05,y05,'linear');
% convert unit from gC m-2 day-1 into gC m-2 yr-1
npp_CMOM05 = npp_CMOM05.*365;
gpp_CMOM05 = gpp_CMOM05.*365;
CUE_CMOM = npp_CMOM05./gpp_CMOM05;
% to calculate the global value
npp_CMOM05PG = npp_CMOM05.* area05.*10^(-15);    % Unit: PgC yr-1
npp_CMOMpg_gb = nansum(npp_CMOM05PG(:));

% merge all observation-based NPP data
NPP_all_obs = [NPP_modisPG_Mgb NPP_CRUpg_Mgb npp_CMOMpg_gb];

% CMIP5 simulations
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

% CMIP6 simulations
cd('H:\CMIP56_Csink\4_MatData\temporal_data\hist_ssp585')
load('2NPP_cmip6_tmp.mat')
NPP6_tmp = [NPPbcc_tmp(2:156),NPPcan_tmp(2:156),NPPcesm_tmp(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   NPPuk_tmp(2:156),NPPipsl_tmp(2:156),NPPmic_tmp(2:156),...          % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   NPPmpi_tmp(2:156),NPPnor_tmp(2:156),...                           % MPI-ESM1-2-LR, NorESM2
   NPPass_tmp(2:156),NPPcnrm_tmp(2:156),NPPec_tmp(2:156)];            % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg
   
NPP5_end5_ag = nanmean(NPP5_tmp(151:155,:),1);
NPP6_end5_ag = nanmean(NPP6_tmp(151:155,:),1); 

NPP_obs_avg = nanmean(NPP_all_obs);

NPP_bias5 = (NPP5_end5_ag - NPP_obs_avg)./NPP_obs_avg;
NPP_bias6 = (NPP6_end5_ag - NPP_obs_avg)./NPP_obs_avg;
NPP_rmse5 = sqrt(mean((NPP5_end5_ag - NPP_obs_avg).^2))
NPP_rmse6 = sqrt(mean((NPP6_end5_ag - NPP_obs_avg).^2))

nanmean(NPP_bias5)
nanmean(NPP_bias6)

clearvars -except cLand_bias5 cLand_bias6 cLand_rmse5 cLand_rmse6 cLand_obs05_Pg...
    NPP_bias5 NPP_bias6 NPP_rmse5 NPP_rmse6 NPP5_end5_ag NPP6_end5_ag NPP_obs_avg...
    CUE_CMOM
%% tuaE = cLand./NPP
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\NOY_circumpolar\NOYpolar_cLand.mat
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\NOY_circumpolar\NOYpolar_GPP56.mat
% global tauE simulated from models
tuaE5_end5_ag = cLand5_M11./NPP5_end5_ag'
tuaE6_end5_ag = cLand6_M11./NPP6_end5_ag'
tuaE_obs = cLand_obs05_Pg./NPP_obs_avg

tuaE_bias5 = (tuaE5_end5_ag - tuaE_obs)./tuaE_obs
tuaE_bias6 = (tuaE6_end5_ag - tuaE_obs)./tuaE_obs
tuaE_rmse5 = sqrt(mean((tuaE5_end5_ag - tuaE_obs).^2))
tuaE_rmse6 = sqrt(mean((tuaE6_end5_ag - tuaE_obs).^2))

nanmean(tuaE_bias5)
nanmean(tuaE_bias6)

clearvars -except cLand_bias5 cLand_bias6 cLand_rmse5 cLand_rmse6 cLand_obs05_Pg...
    NPP_bias5 NPP_bias6 NPP_rmse5 NPP_rmse6 NPP5_end5_ag NPP6_end5_ag NPP_obs_avg...
    CUE_CMOM tuaE_bias5 tuaE_bias6 tuaE_rmse5 tuaE_rmse6

%% Figure3
cmip5_3vars = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','cLand','NPP','tuaE','Nlimitation'});
cmip6_3vars = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','cLand','NPP','tuaE','Nlimitation'});

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
cmip5_3vars.Moldes = Models5';
cmip6_3vars.Moldes = Models6';

cmip5_3vars.cLand = cLand_bias5';
cmip6_3vars.cLand = cLand_bias6';
cmip5_3vars.NPP = NPP_bias5';
cmip6_3vars.NPP = NPP_bias6';
cmip5_3vars.tuaE =  tuaE_bias5;
cmip6_3vars.tuaE =  tuaE_bias6;
NL5 = [0 0 1 0 0 0 0 1 0 0 0]';
NL6 = [0 0 1 1 1 1 1 1  1 0 1]';
cmip5_3vars.Nlimitation = NL5;
cmip6_3vars.Nlimitation = NL6;

[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
figure 
set(gcf,'position',[500,200,885,496])
cmip56(1) = subtightplot(1,2,1,[0.03 0.02],[0.3 0.015],[0.08 0.01])
hold on
set(gca,'Fontname','Arial','FontSize',11)
set(gca,'linewidth',1,'box','off')
%set(gca,'color',[0.92,0.92,0.92])
% CMIP5
h5_X = raincloud_plot(cLand_bias5,'box_on', 0,...
                     'color', [0.97,0.63,0.59], 'cloud_edge_col', [1.00,0.47,0.41],...
                     'bandwidth', 12,'density_type', 'ks','alpha', 0.5);
set(h5_X{2},'Marker','none')                                  
h5_NPP = raincloud_plot(NPP_bias5,'box_on', 0, ...
    'color', [0.80,1.00,0.47], 'cloud_edge_col', [0.56,0.80,0.16],...
    'bandwidth', 12,'density_type', 'ks','alpha', 0.5);
set(h5_NPP{2},'Marker','none')
h5_tuaE = raincloud_plot(tuaE_bias5,'box_on', 0,...
    'color', [0.36,0.92,0.82], 'cloud_edge_col', [0.36,0.92,0.82],...
    'bandwidth', 12,'density_type', 'ks','alpha', 0.2);
set(h5_tuaE{2},'Marker','none')

set(gca, 'YLim',[-8 4],'XLim',[-1.2 1.2]);
plot([0 0], [-8 4],'k--','LineWidth',1.5)
plot( [-1.2 1.2],[0 0],'-','LineWidth',1.8,'color',[0.65 0.65 0.65])


avgcLand5 = nanmean(cmip5_3vars.cLand); 
avgNPP5 = nanmean(cmip5_3vars.NPP); 
avgtuaE5 = nanmean(cmip5_3vars.tuaE); 
for i = 1:11
    
    if cmip5_3vars.Nlimitation(i) == 0 
        leg5_X(i)= plot(cmip5_3vars.cLand(i),-2,'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_NPP(i)= plot(cmip5_3vars.NPP(i),-4.5,'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_tuaE(i)= plot(cmip5_3vars.tuaE(i),-6.5,'Marker','o', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
    
    else 
        leg5_X(i)= plot(cmip5_3vars.cLand(i),-2,'Marker','>', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_NPP(i)= plot(cmip5_3vars.NPP(i),-4.5,'Marker','>', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_tuaE(i)= plot(cmip5_3vars.tuaE(i),-6.5,'Marker','>', 'MarkerEdgeColor', mycolor5(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
    end
    
end

H_pa5_cLand = patch([min(cmip5_3vars.cLand),max(cmip5_3vars.cLand),max(cmip5_3vars.cLand),min(cmip5_3vars.cLand)],[-2.5,-2.5,-1.5,-1.5],cb(4,:));
set(H_pa5_cLand,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
line([avgcLand5 avgcLand5],[-2.5 -1.5],'color','k','linewidth',2);

H_pa5_NPP = patch([min(cmip5_3vars.NPP),max(cmip5_3vars.NPP),max(cmip5_3vars.NPP),min(cmip5_3vars.NPP)],[-5,-5,-4,-4],[0.56,0.80,0.16]);
set(H_pa5_NPP,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
line([avgNPP5 avgNPP5],[-5 -4],'color','k','linewidth',2);

H_pa5_tuaE = patch([min(cmip5_3vars.tuaE),max(cmip5_3vars.tuaE),max(cmip5_3vars.tuaE),min(cmip5_3vars.tuaE)],[-7,-7,-6,-6],[0.00,0.77,0.80]);
set(H_pa5_tuaE,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.35)
line([avgtuaE5 avgtuaE5],[-7 -6],'color','k','linewidth',2);
yticks([-6 -4 -2]);
yticklabels({'\tau_E','NPP','cLand'});
Ylabel = gca
set(Ylabel.YAxis,'FontSize',14,'FontWeight','bold')
set(gca,'Fontname','Arial','FontSize',12)
plot([-1.5, 1.5],[4,4],'k-','LineWidth',1)
text(-1.1,3.3,'(a) CMIP5','Fontname','Arial','FontSize',13)

annotation('arrow',[0.304 0.5289],[0.297,0.297],'HeadStyle','plain',...
    'HeadWidth',10,'HeadLength',9,...
    'LineWidth',4,'Color',[0.97,0.69,0.69])
annotation('arrow',[0.304 0.0805],[0.297,0.297],'HeadStyle','plain',...
    'HeadWidth',10,'HeadLength',9,...
    'LineWidth',4,'Color',[0.48,0.77,0.87])
text(0.7151,-7.53,'Positive bias','Fontname','Arial','FontSize',10,'Color',[0.99,0.03,0.03])
text(-1.1537,-7.5,'Negative bias','Fontname','Arial','FontSize',10,'Color',[0,0,1])

annotation('arrow',[0.2758 0.3038],[0.68 0.68],'HeadStyle','vback1',...
    'HeadWidth',10,'HeadLength',7,...
    'LineWidth',2,'Color',[0,0,0])
text(-0.2681,-0.8853,'Bias','Fontname','Arial','FontSize',12,'Color','r')
annotation('arrow',[0.2758 0.3978],[0.604 0.604],'HeadStyle','vback1',...
    'HeadWidth',10,'HeadLength',7,...
    'LineWidth',2,'Color',[0.89,0.01,0.01])
annotation('arrow',[0.2758 0.2008 ],[0.604 0.604],'HeadStyle','vback1',...
    'HeadWidth',10,'HeadLength',7,...
    'LineWidth',2,'Color',[0.89,0.01,0.01])
text(-0.4175,-3.0029,'Variation','Fontname','Arial','FontSize',12)
tt = text(0,-9.72,'$$\frac{model - Obs}{Obs}$$',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',13)  
set(tt,'Interpreter','latex');
plot([-1.2, -1.2],[-8,4],'k-','LineWidth',1)
plot([1.2, 1.2],[-8,4],'k-','LineWidth',1)

% mark the RMSE for CMIP5
text(-1.16,-1.3,'RMSE:','Fontname','Arial','FontSize',11,'Color','k')
text(-1.14,-2,[num2str(round(cLand_rmse5)),' PgC'],'Fontname','Arial','FontSize',11,'Color','k')
text(-1.14,-4,[num2str(round(NPP_rmse5)),' PgC yr^-^1'],'Fontname','Arial','FontSize',11,'Color','k')
text(-1.14,-6,[num2str(round(tuaE_rmse5)),' yr'],'Fontname','Arial','FontSize',11,'Color','k')
%rectangle('Position',[-1.17,-4,0.2,0.8],'LineWidth',1,'EdgeColor','k');

leg_models5 = legend(leg5_X,Models5,'NumColumns',3)
set(leg_models5,'FontName','Arial','FontSize',8,'Position',[0.069,0.0212,0.42,0.113],...
    'color','none','EdgeColor','k')
text(-1.2518,-10.47,'CMIP5 models:','Fontname','Arial','FontSize',10,'Color','k')


% CMIP6
% 
cmip56(2) = subtightplot(1,2,2,[0.03 0.02],[0.3 0.015],[0.08 0.01])
hold on

set(gca,'Fontname','Arial','FontSize',11)
set(gca,'linewidth',1,'box','off')
set(gca,'Ycolor','k')
h6_X = raincloud_plot(cLand_bias6,'box_on', 0, ...
    'color', [0.97,0.63,0.59], 'cloud_edge_col', [1.00,0.47,0.41],...
    'bandwidth', 12,'density_type', 'ks','alpha', 0.5);
set(h6_X{2},'Marker','none') 
h6_GPP = raincloud_plot(NPP_bias6,'box_on', 0,...
    'color', [0.80,1.00,0.47], 'cloud_edge_col', [0.56,0.80,0.16],...
    'bandwidth', 12,'density_type', 'ks','alpha', 0.5);
set(h6_GPP{2},'Marker','none')
h6_tuaE = raincloud_plot(tuaE_bias6,'box_on', 0,...
    'color', [0.36,0.92,0.82], 'cloud_edge_col', [0.36,0.92,0.82],...
    'bandwidth', 12,'density_type', 'ks','alpha', 0.2);
set(h6_tuaE{2},'Marker','none')

plot([0 0], [-8 4],'k--','LineWidth',1.5)
plot( [-1.5 1.5],[0 0],'-','LineWidth',1.8,'color',[0.65 0.65 0.65])

avgcLand6 = nanmean(cmip6_3vars.cLand); 
avgNPP6 = nanmean(cmip6_3vars.NPP); 
avgtuaE6 = nanmean(cmip6_3vars.tuaE); 
for i = 1:11
    
    if cmip6_3vars.Nlimitation(i) == 0 
        leg6_X(i)= plot(cmip6_3vars.cLand(i),-2,'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg6_NPP(i)= plot(cmip6_3vars.NPP(i),-4,'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg6_tuaE(i)= plot(cmip6_3vars.tuaE(i),-6,'Marker','o', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
    
    else 
        leg6_X(i)= plot(cmip6_3vars.cLand(i),-2,'Marker','>', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg6_NPP(i)= plot(cmip6_3vars.NPP(i),-4,'Marker','>', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg6_tuaE(i)= plot(cmip6_3vars.tuaE(i),-6,'Marker','>', 'MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
    end
    
end

set(gca, 'YLim',[-8 4],'XLim',[-1.2 1.2]);

H_pa5_cLand = patch([min(cmip6_3vars.cLand),max(cmip6_3vars.cLand),max(cmip6_3vars.cLand),min(cmip6_3vars.cLand)],[-2.5,-2.5,-1.5,-1.5],cb(4,:));
set(H_pa5_cLand,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
line([avgcLand6 avgcLand6],[-2.5 -1.5],'color','k','linewidth',2);

H_pa5_NPP = patch([min(cmip6_3vars.NPP),max(cmip6_3vars.NPP),max(cmip6_3vars.NPP),min(cmip6_3vars.NPP)],[-4.5,-4.5,-3.5,-3.5],[0.56,0.80,0.16]);
set(H_pa5_NPP,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
line([avgNPP6 avgNPP6],[-4.5 -3.5],'color','k','linewidth',2);

H_pa5_tuaE = patch([min(cmip6_3vars.tuaE),max(cmip6_3vars.tuaE),max(cmip6_3vars.tuaE),min(cmip6_3vars.tuaE)],[-6.5,-6.5,-5.5,-5.5],[0.00,0.77,0.80]);
set(H_pa5_tuaE,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
line([avgtuaE6 avgtuaE6],[-6.5 -5.5],'color','k','linewidth',2);

yticks([-6 -4 -2]);
yticklabels([]);

plot([-1.2, -1.2],[-8,4],'k-','LineWidth',1)
plot([1.2, 1.2],[-8,4],'k-','LineWidth',1)
plot([-1.2, 1.2],[4,4],'k-','LineWidth',1)
text(-1.1,3.3,'(b) CMIP6','Fontname','Arial','FontSize',13)

annotation('arrow',[0.7667 0.9936],[0.297,0.297],'HeadStyle','plain',...
    'HeadWidth',10,'HeadLength',9,...
    'LineWidth',4,'Color',[0.97,0.69,0.69])
annotation('arrow',[0.7667 0.5452],[0.297,0.297],'HeadStyle','plain',...
    'HeadWidth',10,'HeadLength',9,...
    'LineWidth',4,'Color',[0.48,0.77,0.87])
text(0.7151,-7.53,'Positive bias','Fontname','Arial','FontSize',10,'Color',[0.99,0.03,0.03])
text(-1.1537,-7.5,'Negative bias','Fontname','Arial','FontSize',10,'Color',[0,0,1])
tt = text(0,-9.72,'$$\frac{model - Obs}{Obs}$$',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',13)  
set(tt,'Interpreter','latex');
% mark the RMSE for CMIP6
text(-1.16,-1.3,'RMSE:','Fontname','Arial','FontSize',11,'Color','k')
text(-1.14,-2,[num2str(round(cLand_rmse6)),' PgC'],'Fontname','Arial','FontSize',11,'Color','k')
text(-1.14,-4,[num2str(round(NPP_rmse6)),' PgC yr^-^1'],'Fontname','Arial','FontSize',11,'Color','k')
text(-1.14,-6,[num2str(round(tuaE_rmse6)),' yr'],'Fontname','Arial','FontSize',11,'Color','k')
 

leg_models6 = legend(leg6_X,Models6,'NumColumns',3)
set(leg_models6,'FontName','Arial','FontSize',8,'Position',[0.5411,0.0212,0.4445,0.113],...
    'color','none','EdgeColor','k')
text(-1.2518,-10.47,'CMIP6 models:','Fontname','Arial','FontSize',10,'Color','k')





















