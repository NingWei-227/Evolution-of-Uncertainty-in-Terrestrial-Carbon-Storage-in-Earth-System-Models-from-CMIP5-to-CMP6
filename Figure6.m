clear;clc;
cd('E:\1_Mycase\3_CMIP56_Cland\7_submitNG_sharing\1_Figure3\her_part\CMIP5_temp\matData')

% CMIP5 
load('cLand_cmip5_tmp.mat')
load('NPP_cmip5_tmp.mat')

TauEbcc_tmp =  CLbcc_tmp./NPPbcc_tmp;
TauEbnu_tmp = CLbnu_tmp./NPPbnu_tmp;
TauEcan_tmp = CLcan_tmp./NPPcan_tmp;
TauEccsm_tmp = CLccsm_tmp./NPPccsm_tmp;
TauEgf_tmp = CLgf_tmp./NPPgf_tmp;
TauEhad_tmp = CLhad_tmp./NPPhad_tmp;
TauEipsl_tmp = CLipsl_tmp./NPPipsl_tmp;
TauEmiroc_tmp = CLmiroc_tmp./NPPmiroc_tmp;
TauEmpi_tmp = CLmpi_tmp./NPPmpi_tmp;
TauEmri_tmp = CLmri_tmp./NPPmri_tmp;
TauEnor_tmp = CLnor_tmp./NPPnor_tmp;

save('E:\1_Mycase\3_CMIP56_Cland\8_CodeFigures_JClimate\MainFigures\M_Figure6\Figure6_version2\her_part\CMIP5_temp\matData\tuaE_cmip5_tmp.mat',...
    'TauEbcc_tmp','TauEbnu_tmp','TauEcan_tmp','TauEccsm_tmp','TauEgf_tmp','TauEhad_tmp','TauEipsl_tmp','TauEmiroc_tmp','TauEmpi_tmp','TauEmri_tmp','TauEnor_tmp');

clear;clc;
% CMIP6
cd('E:\1_Mycase\3_CMIP56_Cland\8_CodeFigures_JClimate\MainFigures\M_Figure6\Figure6_version2\her_part\CMIP6_temp\matData')
load('2cLand_cmip6_tmp.mat')
load('2NPP_cmip6_tmp.mat')

TauEass_tmp = CLass_tmp./NPPass_tmp;
TauEbcc_tmp = CLbcc_tmp./NPPbcc_tmp;
TauEcan_tmp = CLcan_tmp./NPPcan_tmp;

TauEcesm1m_tmp = CLcesm1m_tmp./NPPcesm_tmp(1:165);
TauEcesm_tmp = CLcesm_tmp./NPPcesm_tmp;
TauEcesm_W_tmp = CLcesm_W_tmp./NPPcesm_W_tmp;

TauEcnrm_tmp = CLcnrm_tmp./NPPcnrm_tmp;
TauEec_tmp = CLec_tmp./NPPec_tmp;
TauEipsl_tmp = CLipsl_tmp./NPPipsl_tmp;

TauEmic_tmp = CLmic_tmp./NPPmic_tmp;
TauEmpi_tmp = CLmpi_tmp./NPPmpi_tmp;
TauEnor1m_tmp = CLnor1m_tmp./NPPnor_tmp(1:165);

TauEnor_tmp = CLnor_tmp./NPPnor_tmp;
TauEuk_tmp = CLuk_tmp./NPPuk_tmp;

save('E:\1_Mycase\3_CMIP56_Cland\8_CodeFigures_JClimate\MainFigures\M_Figure6\Figure6_version2\her_part\CMIP6_temp\matData\2tuaE_cmip6_tmp.mat',...
    'TauEass_tmp','TauEbcc_tmp','TauEcan_tmp','TauEcesm1m_tmp','TauEcesm_tmp','TauEcesm_W_tmp','TauEcnrm_tmp','TauEec_tmp','TauEipsl_tmp',...
    'TauEmic_tmp','TauEmpi_tmp','TauEnor1m_tmp','TauEnor_tmp','TauEuk_tmp')








