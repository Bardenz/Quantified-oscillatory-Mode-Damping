%% EMT Simu (get Measurents for SS values Z responses) - Analytical Modelling - Validation (Using EMT Measurents)
addpath([pwd '/lib_para'])

my_MMC_param

global SEQUENCE_NO ANALYTICAL_FLAG MEAS_FLAG
SEQUENCE_NO = 1; ANALYTICAL_FLAG = 1; MEAS_FLAG = 1;

%% General MMC settings
MMC.CtrlMode = GFL_APC; %% selection of control mode: GFM_VF = 0; GFL_APC = 1; GFL_DVC = 2; GFM_APC = 3; (values programmed in EMT models)
MMC.Ctrl.Modulation = 1; % 1 UdcN OR Udc0, 2 udc (or averaged udc??), 3, arm capacitor voltage sum
MMC.Ctrl.CCSC.Flag = 1; %% 0 for disabled, 1 for closed-loop control, 2 for open-loop control - not used anymore since 2024
MMC.Ctrl.ZSCC.Flag = 0; %% 0 for disabled, 1 for energy control, 2 for Idc active damping (active if CCSC_FLAG = 1)
MMC.Ctrl.IsPacFFinWt = 1; %% 1 denotes having power feedforward in energy control, 0 for not
%% Settings for GFL
MMC.Ctrl.IsCrossGFL = 0; %% 1 for cross control, 0 for classic control
MMC.Ctrl.IsVacGFL = 0; %% 1 for PV GFL mode, 0 for PQ GFL mode
%% Settings for GFM
MMC.Ctrl.IsCCinGFM = 0; %% 0 for without inner loop, 1 for with inner loop
MMC.Ctrl.IsVIinGFM = 0; %% 1 for having virtual impedance (VI), and 0 for not
MMC.Ctrl.IsPfDroopGFM = 0; %% as default only enabled beyond [49.8 50.2] Hz in EMT (to check whether it is commented through)
%% Setpoints
MMC.Set.Pref = 0.6; %1001e6/VSC1.pu.Sb;
MMC.Set.Qref = 0.0; %0*61e6/VSC1.pu.Sb;

%% create MMC ojects
MMC1 = MMC;
MMC1.CtrlMode = GFM_APC;
MMC1.Ctrl.IsCCinGFM = 1;
MMC1.Ctrl.IsPfDroopGFM = 0;
MMC1.Ctrl.IsVIinGFM = 0; 
%
MMC2 = MMC;
MMC2.CtrlMode = GFL_DVC;
%MMC2.ZSCC.Flag = 1;
%
MMC3 = MMC;
MMC3.CtrlMode = GFM_APC;
MMC3.Set.Pref = 0.6;