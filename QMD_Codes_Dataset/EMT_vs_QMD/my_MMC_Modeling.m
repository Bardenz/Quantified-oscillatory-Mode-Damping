%
%% EMT Simu (get Measurents for SS values Z responses) - Analytical Modelling - Validation (Using EMT Measurents)
%
function YmmcStru = my_MMC_Modeling(mmcSetup)

global SEQUENCE_NO ANALYTICAL_FLAG MEAS_FLAG

addpath([pwd '/lib'])

%GridParallel_CC_init_Inverter

my_MMC_param

currentFolder = pwd;
SUB_FOLDER_STR = '\MeasData\m_Ypndc_';
Z_FOLDER_STR = '\Zf_Data\m_Ypndc_';
kHz_val = sprintf('%0.0f',1/Ts_Power/1000);

%% see my_MMC_param.m for control mode definition
if mmcSetup.CtrlMode == GFM_VF
    measFolder = strcat(currentFolder,SUB_FOLDER_STR,'GFM_VF_',kHz_val,'kHz');
    zFolder = strcat(currentFolder,Z_FOLDER_STR,'GFM_VF_',kHz_val,'kHz');
    FILE_PREFIX = '\GFM_VF_';
elseif mmcSetup.CtrlMode == GFL_APC
    measFolder = strcat(currentFolder,SUB_FOLDER_STR,'GFL_APC_',kHz_val,'kHz');
    zFolder = strcat(currentFolder,Z_FOLDER_STR,'GFL_APC_',kHz_val,'kHz');
    FILE_PREFIX = '\GFL_APC_'; %% adding woVff_ if voltage feedforward is removed
elseif mmcSetup.CtrlMode == GFL_DVC
    measFolder = strcat(currentFolder,SUB_FOLDER_STR,'GFL_DVC_',kHz_val,'kHz');
    zFolder = strcat(currentFolder,Z_FOLDER_STR,'GFL_DVC_',kHz_val,'kHz');
    FILE_PREFIX = '\GFL_DVC_';
elseif mmcSetup.CtrlMode == GFM_APC
    measFolder = strcat(currentFolder,SUB_FOLDER_STR,'GFM_APC_',kHz_val,'kHz');
    zFolder = strcat(currentFolder,Z_FOLDER_STR,'GFM_APC_',kHz_val,'kHz');
    FILE_PREFIX = '\GFM_APC_';
end

Pset= mmcSetup.Set.Pref;
% if mmcSetup.CtrlMode == GFL_DVC
%     Pset= -mmcSetup.Set.Pref;
% else
%     Pset= mmcSetup.Set.Pref;
% end
FILE_PREFIX = strcat(FILE_PREFIX,'P',strrep(num2str(Pset), '.', ''),'pu_');

if(~isfolder(measFolder))
    mkdir(measFolder);
end
if(~isfolder(zFolder))
    mkdir(zFolder);
end

if mmcSetup.Ctrl.Modulation== 2
    FILE_PREFIX = strcat(FILE_PREFIX,'UCM2_');
elseif mmcSetup.Ctrl.Modulation == 3
    FILE_PREFIX = strcat(FILE_PREFIX,'CM3_');
end

if mmcSetup.Ctrl.IsVacGFL == 1 && mmcSetup.CtrlMode ~= GFM_APC
    FILE_PREFIX = strcat(FILE_PREFIX,'V_');
end
if mmcSetup.Ctrl.IsCrossGFL == 1 && mmcSetup.CtrlMode ~= GFM_APC
    FILE_PREFIX = strcat(FILE_PREFIX,'CROSS_');
end
if mmcSetup.Ctrl.IsPacFFinWt == 0 && mmcSetup.CtrlMode == GFM_APC
    FILE_PREFIX = strcat(FILE_PREFIX,'DECOUPLED_');
end
if mmcSetup.Ctrl.IsCCinGFM == 1 && mmcSetup.CtrlMode == GFM_APC
    FILE_PREFIX = strcat(FILE_PREFIX,'CC_');
end
if mmcSetup.Ctrl.IsVIinGFM == 1 && mmcSetup.CtrlMode == GFM_APC
    FILE_PREFIX = strcat(FILE_PREFIX,'VI_');
end

%% GRID STRENGTH TEST
if Psc ~= Pnom*20
    %FILE_PREFIX = strcat(FILE_PREFIX,'SCR', num2str(Psc/Pnom),'_');
    if Psc/Pnom >3 && mmcSetup.CtrlMode == GFM_APC && mmcSetup.Ctrl.IsCCinGFM == 1
        FILE_PREFIX = strcat(FILE_PREFIX,'SCR', num2str(3),'_');
    else
        FILE_PREFIX = strcat(FILE_PREFIX,'SCR', num2str(Psc/Pnom),'_');
    end  
end

    %% first number after CCSC denotes CCSC type, while the second number denotes ZSCC type
    MEAS_FILE_P = strcat(measFolder,FILE_PREFIX,'CCSC',num2str(mmcSetup.Ctrl.CCSC.Flag),num2str(mmcSetup.Ctrl.ZSCC.Flag),'_p.mat');%_add_mdc.mat');
    MEAS_FILE_N = strcat(measFolder,FILE_PREFIX,'CCSC',num2str(mmcSetup.Ctrl.CCSC.Flag),num2str(mmcSetup.Ctrl.ZSCC.Flag),'_n.mat');%_add_mdc.mat');
    MEAS_FILE_DC = strcat(measFolder,FILE_PREFIX,'CCSC',num2str(mmcSetup.Ctrl.CCSC.Flag),num2str(mmcSetup.Ctrl.ZSCC.Flag),'_dc.mat');%_add_mdc.mat');

    IMPEDANCE_FILE = strcat(zFolder,FILE_PREFIX,'CCSC',num2str(mmcSetup.Ctrl.CCSC.Flag),num2str(mmcSetup.Ctrl.ZSCC.Flag),'.mat');%_add_mdc.mat');
    
    %% if impedance file exists, read the data and return
    if isfile(IMPEDANCE_FILE)
        YmmcStru = load(IMPEDANCE_FILE).YmmcStru;
        return;
    end

    %% Run simulations to generate data files
    SEQUENCE_NO = 1; %% default value for p-sequence
    %% return,
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SIM_LOG_FILE=strcat(currentFolder,'\DATA_LOG_FILE.mat');
    if isfile(SIM_LOG_FILE) % returns 1 if fileName is a file located on the specified path
        delete(SIM_LOG_FILE);
    end

    if ~isfile(MEAS_FILE_P) 
        SEQUENCE_NO = 1,
        %% warning('off','all'); %% does not work
        sim('MMC_Modeling_EMT')        
        movefile(SIM_LOG_FILE,MEAS_FILE_P);
    end
    if ~isfile(MEAS_FILE_N)
        SEQUENCE_NO = 2,
        sim('MMC_Modeling_EMT')
        movefile(SIM_LOG_FILE,MEAS_FILE_N);
    end
    %% extra simulations to run for dc perturbation
    if  ~isfile(MEAS_FILE_DC) %&& MMC_Y_PORT_NUM == 3
        SEQUENCE_NO = 0,
        sim('MMC_Modeling_EMT')
        movefile(SIM_LOG_FILE,MEAS_FILE_DC);
    end


%% load data files
Mp = load(MEAS_FILE_P);

MEAS_FILE_P,
Mn = load(MEAS_FILE_N);
MdcDATA = load(MEAS_FILE_DC);
RowColNum = size(Mp.logsout.get);
MeasNum = RowColNum(1);

%THETA_G_EXISTS = 0;
%% GET measured data from logged files
for i=1:MeasNum
    if strcmp(Mp.logsout.get(i).Name,'<meas_u_abc>') == 1
        Meas.Uac_inj1 = Mp.logsout.get(i).Values.Data/Ubase;
        Meas.Uac_inj2 = Mn.logsout.get(i).Values.Data/Ubase;
        Meas.Uac_inj3 = MdcDATA.logsout.get(i).Values.Data/Ubase;
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_i_abc>') == 1
        Meas.Iac_inj1 = -Mp.logsout.get(i).Values.Data/Ibase; % total output current
        Meas.Iac_inj2 = -Mn.logsout.get(i).Values.Data/Ibase;
        Meas.Iac_inj3 = -MdcDATA.logsout.get(i).Values.Data/Ibase;
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_u_dc>') == 1
        Meas.Udc_inj1 = Mp.logsout.get(i).Values.Data/Ubase;
        Meas.Udc_inj2 = Mn.logsout.get(i).Values.Data/Ubase;
        Meas.Udc_inj3 = MdcDATA.logsout.get(i).Values.Data/Ubase;
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_i_dc>') == 1
        Meas.Idc_inj1 = -Mp.logsout.get(i).Values.Data/Ibase; % representation of current direction to be checked
        Meas.Idc_inj2 = -Mn.logsout.get(i).Values.Data/Ibase;
        Meas.Idc_inj3 = -MdcDATA.logsout.get(i).Values.Data/Ibase;    
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_uc_sum>') == 1
        Meas.uc_sum = Mp.logsout.get(i).Values.Data;
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_uc_diff>') == 1
        Meas.uc_diff = Mp.logsout.get(i).Values.Data;
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_icir_abc>') == 1
        Meas.icir_abc = Mp.logsout.get(i).Values.Data;    
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_uac_dq>') == 1
        Meas.uac_dq = Mp.logsout.get(i).Values.Data;
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_iac_dq>') == 1
        Meas.iac_dq = Mp.logsout.get(i).Values.Data;
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_uac_ref_dq>') == 1
        Meas.uac_ref_dq = Mp.logsout.get(i).Values.Data;
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_ucir_ref_dq>') == 1
        Meas.ucir_ref_dq = Mp.logsout.get(i).Values.Data;
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_theta_ref>') == 1
        Meas.theta = Mp.logsout.get(i).Values.Data;
    elseif strcmp(Mp.logsout.get(i).Name,'<meas_theta_GFM>') == 1
        Meas.theta_GFM = Mp.logsout.get(i).Values.Data;
    % elseif strcmp(Mp.logsout.get(i).Name,'<meas_theta_g>') == 1
    %     Meas.theta_g = Mp.logsout.get(i).Values.Data;
    %     %THETA_G_EXISTS = 1;
    end
end
%% END
YmmcStru.SimuMeas = Meas;

if MEAS_FLAG == 0
    mMeas = [];
    YmmcStru.mMeas = [];
else
    mMeas = my_MEAS_fun(Meas);
    YmmcStru.mMeas = mMeas;
end

if ANALYTICAL_FLAG == 0
    YmmcStru.f_buf = [];
    YmmcStru.Zmmc_pu_buf = [];
    YmmcStru.Ymmc_pu_buf = [];
    return;
end
%% Getting Steady-State Values (Here control and power sampling frequencies should be distinguished)
if mmcSetup.CtrlMode == GFM_VF
    theta0 = 0;
elseif mmcSetup.CtrlMode == GFM_APC
    theta0 = Meas.theta_GFM(round((20-1)/Ts_Control));
    theta0_PLL = Meas.theta(round((20-1)/Ts_Control));
else
    theta0 = Meas.theta(round((20-1)/Ts_Control));
end
% if length(Meas.mdc_dq) >= round((20-1)/Ts_Control)
%     m0_cir=Meas.ucir_ref_dq(round((20-1)/Ts_Control),:).';
% else
%     m0_cir=Meas.mdc_dq(1,:).';
% end
Uref_ac=Meas.uac_ref_dq(round((20-1)/Ts_Control),:).';
Uref_cir=Meas.ucir_ref_dq(round((20-1)/Ts_Control),:).';
%% from AC/DC Uref to mac/mdc considering modulation compensation, but NOT INCLUDING the transformation from ac/dc to u/l
m0_ac = Uref_ac;
m0_cir = Uref_cir;

theta0 = theta0 - PhaseDelay*0;

u0_ac=Meas.uac_dq(round((20-1)/Ts_Control),:).';
i0_ac=Meas.iac_dq(round((20-1)/Ts_Control),:).';

U1d=u0_ac(1); U1q=u0_ac(2); %%% Us_vec=U1d+jU1q, measured values
I1d=i0_ac(1); I1q=i0_ac(2); %%% Ig1_vec=Ig1d+jIg1q, measured values
M1d=m0_ac(1); M1q=m0_ac(2);

%% following fft results have been checked with Simulink Sequence-Analyzer block (angle of u0_abc also compared with PLL angle)
Uac_fft = my_fft_fun(Meas.Uac_inj1,Ts_Power);
Iac_fft = my_fft_fun(-Meas.Iac_inj1,Ts_Power);
Icir_fft = my_fft_fun(Meas.icir_abc,Ts_Power); %% p.u. values are measured
Uc_sum_fft = my_fft_fun(Meas.uc_sum,Ts_Power);
Uc_diff_fft = my_fft_fun(Meas.uc_diff,Ts_Power);
% mag=abs(Icir2_abc),
% deg=angle(Icir2_abc)*180/pi,

YmmcStru.m0_ac = m0_ac;
YmmcStru.m0_cir = m0_cir;

%% Analytical Impedance modelling - since March 2023

%% Power module model derived using HSS method by setting the considered harmonic order as h=3
%% Case 1 - open loop: setting m_cir=0, and m_ac is fixed for generating a constant V(f).
Tabc2pn_hss = blkdiag(Tabc2pn, Tabc2pn, Tabc2pn, Tabc2pn, Tabc2pn, Tabc2pn, Tabc2pn);
Tpn2abc_hss = blkdiag(Tpn2abc, Tpn2abc, Tpn2abc, Tpn2abc, Tpn2abc, Tpn2abc, Tpn2abc);

Tabc2pnz_hss = blkdiag(Tabc2pnz, Tabc2pnz, Tabc2pnz, Tabc2pnz, Tabc2pnz, Tabc2pnz, Tabc2pnz);
Tpnz2abc_hss = blkdiag(Tpnz2abc, Tpnz2abc, Tpnz2abc, Tpnz2abc, Tpnz2abc, Tpnz2abc, Tpnz2abc);

%Tpn2dq_hss = blkdiag(Tpn2dq, Tpn2dq, Tpn2dq, Tpn2dq, Tpn2dq, Tpn2dq, Tpn2dq);
F_shift = diag([-3j*ones(1,12)	-2j*ones(1,12)	-1j*ones(1,12)	0j*ones(1,12)   1j*ones(1,12)   2j*ones(1,12)   3j*ones(1,12)])*wb;

%% Mdc0 is 1 at least for Modulation mode 1 (Udc) and 2 (udc), suitability for Compensated Modulation to check
Mdc0 = ones(3,1); %% Udc0_ref=Udc0 => mdc0=Udc0_ref/Udc0=1;

Mac1_dq = [m0_ac(1)    m0_ac(2)].';

Ttheta0 = [cos(theta0)      sin(theta0)
          -sin(theta0)      cos(theta0)];        
Ttheta_0 = [cos(-theta0)        sin(-theta0)
            -sin(-theta0)       cos(-theta0)];    
%% following representations have been checked with mac_abc measurements 
Mac1_pn = inv(Tpn2dq)*Ttheta_0*Mac1_dq; %% positive sequence complex phasors
%% Be careful with the transformation between trigonometric and exponential representations
Mac1 = Tpn2abc*[Mac1_pn(1) 0].'; %% abc complex phasors after spliting the +/- f1 frequency Fourier coefficients
Mac_1 = Tpn2abc*[0 Mac1_pn(2)].'; %% %% abc complex phasors

Mdc2_dq = [m0_cir(1)     m0_cir(2)].';
Mdc2_theta0 = -(-2*theta0)+pi/2*0;
T_2theta0 = [cos(Mdc2_theta0)      sin(Mdc2_theta0)
            -sin(Mdc2_theta0)      cos(Mdc2_theta0)];
Mdc2_pn = inv(Tpn2dq)*T_2theta0*Mdc2_dq; %% including the devision of 2 for getting complex fourier coefficients
%% Be careful with the transformation between trigonometric and exponential representations
Mdc2 = 1*Tpn2abc*[Mdc2_pn(1)    0].'; %% abc
Mdc_2 = 1*Tpn2abc*[0    Mdc2_pn(2)].'; %% abc
%% Figure out Why m2 insersion indexes strongly vary the three-phase balance of iac steady-state values 
if mmcSetup.Ctrl.CCSC.Flag == 0
    Mdc2 = zeros(3,1);
    Mdc_2 = zeros(3,1);
end
    
if IS_VNO_INCLUDED == 1
    Mdc0A = -[-1/3 1/6 1/6; 1/6 -1/3 1/6; 1/6 1/6 -1/3]*diag(Mdc0)*2;
    Mac1A = [2/3 -1/3 -1/3; -1/3 2/3 -1/3; -1/3 -1/3 2/3]*diag(Mac1);
    Mac_1A = [2/3 -1/3 -1/3; -1/3 2/3 -1/3; -1/3 -1/3 2/3]*diag(Mac_1);
    Mdc2A = -[-1/3 1/6 1/6; 1/6 -1/3 1/6; 1/6 1/6 -1/3]*diag(Mdc2)*2;
    Mdc_2A = -[-1/3 1/6 1/6; 1/6 -1/3 1/6; 1/6 1/6 -1/3]*diag(Mdc_2)*2;
else
    Mdc0A = diag(Mdc0);
    Mac1A = diag(Mac1);
    Mac_1A = diag(Mac_1);
    Mdc2A = diag(Mdc2);
    Mdc_2A = diag(Mdc_2);
end
%% A_h is a 12-by-12 coefficient matrix for x = [iac_abc   icir*ones(1,3)    uc_abc_cm   uc_abc_dm]
%% Figure out why Mdc0 should be set as 1 while there is no such component in the control when using the default PWM modulation
A0 = [  -Rarm_pu/(Larm_pu/wb)*eye(3)        zeros(3)                                zeros(3)                        -Mdc0A/2/(Larm_pu/wb)
        zeros(3)                            -Rarm_pu/(Larm_pu/wb)*eye(3)            -diag(Mdc0)/4/(Larm_pu/wb)      zeros(3)
        zeros(3)                            Nb_PM/(Csm_pu/wb)*diag(Mdc0)            zeros(3)                        zeros(3)
        Nb_PM/2/(Csm_pu/wb)*diag(Mdc0)      zeros(3)                                zeros(3)                        zeros(3)];

A1 = [  zeros(3)                            zeros(3)                                Mac1A/(Larm_pu/wb)              zeros(3)
        zeros(3)                            zeros(3)                                zeros(3)                        diag(Mac1)/2/(Larm_pu/wb)
        -Nb_PM/(Csm_pu/wb)*diag(Mac1)       zeros(3)                                zeros(3)                        zeros(3)
        zeros(3)                            -2*Nb_PM/(Csm_pu/wb)*diag(Mac1)         zeros(3)                        zeros(3)];
    
A_1 = [ zeros(3)                            zeros(3)                                Mac_1A/(Larm_pu/wb)             zeros(3)
        zeros(3)                            zeros(3)                                zeros(3)                        diag(Mac_1)/2/(Larm_pu/wb)
        -Nb_PM/(Csm_pu/wb)*diag(Mac_1)      zeros(3)                                zeros(3)                        zeros(3)
        zeros(3)                            -2*Nb_PM/(Csm_pu/wb)*diag(Mac_1)        zeros(3)                        zeros(3)];

A2 = [  zeros(3)                            zeros(3)                                zeros(3)                        -Mdc2A/2/(Larm_pu/wb)
        zeros(3)                            zeros(3)                                -Mdc2/4/(Larm_pu/wb).*eye(3)	zeros(3)
        zeros(3)                            Nb_PM/(Csm_pu/wb)*Mdc2.*eye(3)          zeros(3)                        zeros(3)
        Nb_PM/2/(Csm_pu/wb)*Mdc2.*eye(3)	zeros(3)                                zeros(3)                        zeros(3)];
    
A_2 = [ zeros(3)                            zeros(3)                                zeros(3)                        -Mdc_2A/2/(Larm_pu/wb)
        zeros(3)                            zeros(3)                                -Mdc_2/4/(Larm_pu/wb).*eye(3)	zeros(3)
        zeros(3)                            Nb_PM/(Csm_pu/wb)*Mdc_2.*eye(3)         zeros(3)                        zeros(3)
        Nb_PM/2/(Csm_pu/wb)*Mdc_2.*eye(3)	zeros(3)                                zeros(3)                        zeros(3)];

%A2 = zeros(12); A_2 = zeros(12);
A3 = zeros(12); A_3 = zeros(12);

%% A_hss is a 12x7-by-12x7 matrix for steady-state coefficients of abc-variables
%% x = [iac_abc   icir*ones(1,3)    uc_abc_cm   uc_abc_dm]
A_Tpz =[A0          A_1         A_2         A_3         zeros(12)       zeros(12)       zeros(12)
        A1          A0          A_1         A_2         A_3             zeros(12)       zeros(12)
        A2          A1          A0          A_1         A_2             A_3             zeros(12)
        A3          A2          A1          A0          A_1             A_2             A_3
        zeros(12)   A3          A2          A1          A0              A_1             A_2
        zeros(12)	zeros(12)	A3          A2          A1              A0              A_1
        zeros(12)	zeros(12)	zeros(12)	A3          A2              A1              A0  ];
%% In case of u = [uac_abc     udc*ones(1,3)].', Bh/B_h = zeros(12,6)
%% Note that the power-stage model is formulated under per-unit system

B0 = [-2/(Larm_pu/wb).*eye(3)	zeros(3)
      zeros(3)                  1/2/(Larm_pu/wb).*eye(3)
      zeros(3)                  zeros(3)
      zeros(3)                  zeros(3)    ];
B1 = zeros(12,6); B2 = zeros(12,6); B3 = zeros(12,6); B_1 = zeros(12,6); B_2 = zeros(12,6); B_3 = zeros(12,6);
%% B_hss is a (7*12)x(7*6) matrix
B_Tpz =[B0              B_1             B_2             B_3         zeros(12,6)     zeros(12,6)     zeros(12,6)
        B1              B0              B_1             B_2         B_3             zeros(12,6)     zeros(12,6)
        B2              B1              B0              B_1         B_2             B_3             zeros(12,6)
        B3              B2              B1              B0          B_1             B_2             B_3
        zeros(12,6)     B3              B2              B1          B0              B_1             B_2
        zeros(12,6)     zeros(12,6)     B3              B2          B1              B0              B_1
        zeros(12,6)     zeros(12,6)     zeros(12,6)     B3          B2              B1              B0  ];

%% -------------------------------------------------------------------------------------------
%% Steady-State values:Udc0 calculated using Ubase instead of Udcb, which however does not 
u_SS0 = [zeros(3,1); Vnom_dc/Ubase*ones(3,1)]; 
%% Note that in p.u. representation the steady-state udc is equal to 2*Udc0 when applying Ubase_dc = 2*Ubase_ac as dc base
u_SS1 = [0.5*Uac_fft.f1; zeros(3,1)]; 
u_SS_1 = [0.5*conj(Uac_fft.f1); zeros(3,1)];
u_SSh = [zeros(6,1); zeros(6,1); u_SS_1; u_SS0; u_SS1; zeros(6,1); zeros(6,1)];
if 0    
%% solving the multi-frequency state-space equation at steady state =>> does not give correct results yet
x_SS = -inv(A_Tpz - 1*F_shift) * (B_Tpz*u_SSh);
else
%% ------------------- estimated and measured values ---------------------
%% x = [iac_abc   icir*ones(1,3)    uc_abc_sum   uc_abc_diff].'
x_SS0 = [Iac_fft.f0; Icir_fft.f0; Uc_sum_fft.f0; Uc_diff_fft.f0];
x_SS1 = 0.5*[Iac_fft.f1; Icir_fft.f1; Uc_sum_fft.f1; Uc_diff_fft.f1];
x_SS_1 = conj(x_SS1);
x_SS2 = 0.5*[Iac_fft.f2; Icir_fft.f2; Uc_sum_fft.f2; Uc_diff_fft.f2];
x_SS_2 = conj(x_SS2);
x_SS3 = 0.5*[Iac_fft.f3; Icir_fft.f3; Uc_sum_fft.f3; Uc_diff_fft.f3];
x_SS_3 = conj(x_SS3);
x_SS = [x_SS_3; x_SS_2; x_SS_1; x_SS0; x_SS1; x_SS2; x_SS3];
end
%% Call the function of C_SS_fun()
Css = my_MMC_Css_fun(x_SS);
%% 12x7 by 6x7, coefficient of [mac     mdc]
C_Tpz =[Css.C0          Css.C_1         Css.C_2         Css.C_3         zeros(12,6)     zeros(12,6)     zeros(12,6)
        Css.C1          Css.C0          Css.C_1         Css.C_2         Css.C_3         zeros(12,6)     zeros(12,6)
        Css.C2          Css.C1          Css.C0          Css.C_1         Css.C_2         Css.C_3         zeros(12,6)
        Css.C3          Css.C2          Css.C1          Css.C0          Css.C_1         Css.C_2         Css.C_3
        zeros(12,6)     Css.C3          Css.C2          Css.C1          Css.C0          Css.C_1         Css.C_2
        zeros(12,6)     zeros(12,6)     Css.C3          Css.C2          Css.C1          Css.C0          Css.C_1
        zeros(12,6)     zeros(12,6)     zeros(12,6)     Css.C3          Css.C2          Css.C1          Css.C0  ];
%% x_abc = [iac_abc   icir*ones(1,3)    uc_abc_cm   uc_abc_dm].' 
%% u_abc = [uac_abc     udc*ones(1,3)].'
IN.u0_ac = u0_ac/1;
IN.i0_ac = i0_ac/1;
IN.m0_ac = Uref_ac;
IN.theta0 = theta0/1;
IN.m0_cir = Uref_cir;
IN.i0_cir2 = T_2theta0*Tpn2dq*Tabc2pn*Icir_fft.f2;
%% the dq steady-state values below still need to be checked
%% the real()function happens to capature controller dq-frame steady-state values, which exihbits negligible impedance-shaping effects
IN.uc_sum0_f0z = Uc_sum_fft.f0(1);
IN.uc_sum0_f2dq = inv(T_2theta0)*Tpn2dq*Tabc2pn*Uc_sum_fft.f2;
IN.uc_diff0_f1dq = Ttheta0*Tpn2dq*Tabc2pn*Uc_diff_fft.f1;
%     IN.uc_sum0_f2dq = real(inv(T_2theta0)*Tpn2dq*Tabc2pn*Uc_sum_fft.f2);
%     IN.uc_diff0_f1dq = real(Ttheta0*Tpn2dq*Tabc2pn*Uc_diff_fft.f1);

Tabc2pnz_ss = blkdiag(Tabc2pnz, Tabc2pnz, Tabc2pnz, Tabc2pnz);
Tabc2pnz_hss_x = blkdiag(Tabc2pnz_ss, Tabc2pnz_ss, Tabc2pnz_ss, Tabc2pnz_ss, Tabc2pnz_ss, Tabc2pnz_ss, Tabc2pnz_ss);
Tpnz2abc_ss = blkdiag(Tpnz2abc, Tpnz2abc);
Tpnz2abc_hss_u = blkdiag(Tpnz2abc_ss, Tpnz2abc_ss, Tpnz2abc_ss, Tpnz2abc_ss, Tpnz2abc_ss, Tpnz2abc_ss, Tpnz2abc_ss);

%% GET analytical Impedance/Admittance BY frequency-scanning
Cnt = 0;
for k=Hmin/Hstep:1:Hmax/Hstep
%for k=1-1/Hstep:1:Hmax/Hstep+1/Hstep
    w=k*Hstep*wb;
    s=j*w;
    %% h is exclusively used for denoting harmonic orders of MMC model

    if w >= 2*pi*1 && w < 2*pi*10 && rem(w,2*pi*0.01)~=0
        continue;
    elseif w >= 2*pi*10 && w < 2*pi*100 && (rem(w,2*pi*0.1)~=0 || w==wb) %% Remove f1 component to prevent ERROR IN Zm frequency sweep
        continue;
    elseif w >= 2*pi*100 && (rem(w,2*pi*1)~=0 || w==2*wb || w==3*wb) %% Remove f1 component to prevent ERROR IN Zm frequency sweep
        continue;
    end

    Cnt = Cnt + 1;
    f_buf(Cnt) = w/2/pi;
    
    C2x2=[0      1; -1     0];
    %% Refer to Chen's paper "MMC impedance modelling ... in close proximity"  
    %% 6x7 by 12x7
    IN.s = s;
    GmxStru = my_MMC_Gmx_fun(IN,mmcSetup);
    %Gmxh = blkdiag(Gmx_fun(s,-3), Gmx_fun(s,-2), Gmx_fun(s,-1), Gmx_fun(s,0), Gmx_fun(s,1), Gmx_fun(s,2), Gmx_fun(s,3));
    Gmxh = [GmxStru.f_3.D	    zeros(6,12)         GmxStru.f_1.R	        zeros(6,12)         zeros(6,12)             zeros(6,12)         zeros(6,12)
            zeros(6,12)         GmxStru.f_2.D	    zeros(6,12)             GmxStru.f0.R        GmxStru.f1.ccsc_3*1     zeros(6,12)         zeros(6,12)
            GmxStru.f_3.L	    zeros(6,12)         GmxStru.f_1.D	        zeros(6,12)         GmxStru.f1.R            zeros(6,12)         zeros(6,12)
            zeros(6,12)         GmxStru.f_2.L	    zeros(6,12)             GmxStru.f0.D        GmxStru.f1.W1*0	        GmxStru.f2.R        zeros(6,12)
            zeros(6,12)         zeros(6,12)         GmxStru.f_1.L	        zeros(6,12)         GmxStru.f1.D            zeros(6,12)         GmxStru.f3.R
            zeros(6,12)         zeros(6,12)         GmxStru.f_1.ccsc3*1     GmxStru.f0.L        zeros(6,12)             GmxStru.f2.D        zeros(6,12)
            zeros(6,12)         zeros(6,12)         zeros(6,12)             zeros(6,12)         GmxStru.f1.L            zeros(6,12)         GmxStru.f3.D ];
    
    GmuStru = my_MMC_Gmu_fun(IN,mmcSetup);
    %% Refer to Chen's paper "MMC impedance modelling ... in close proximity" 
    if 1
    Gmuh = [GmuStru.f_3.D	    zeros(6)            GmuStru.f_1.R	        zeros(6)            zeros(6)                zeros(6)        zeros(6)
            zeros(6)            GmuStru.f_2.D	    zeros(6)                GmuStru.f0.R        GmuStru.f1.ccsc_3*1     zeros(6)        zeros(6)
            GmuStru.f_3.L	    zeros(6)            GmuStru.f_1.D	        GmuStru.f0.dc_1	    GmuStru.f1.R            zeros(6)        zeros(6)
            zeros(6)            GmuStru.f_2.L	    zeros(6)                GmuStru.f0.D        zeros(6)                GmuStru.f2.R	zeros(6)
            zeros(6)            zeros(6)            GmuStru.f_1.L	        GmuStru.f0.dc1	    GmuStru.f1.D            zeros(6)        GmuStru.f3.R
            zeros(6)            zeros(6)            GmuStru.f_1.ccsc3*1     GmuStru.f0.L        zeros(6)                GmuStru.f2.D	zeros(6)
            zeros(6)            zeros(6)            zeros(6)                zeros(6)            GmuStru.f1.L            zeros(6)        GmuStru.f3.D ];
    else
    Gmuh = [GmuStru.f_3.D	    GmuStru.f_2.dc_1	GmuStru.f_1.R	    zeros(6)            zeros(6)            zeros(6)            zeros(6)
            GmuStru.f_3.dc1	    GmuStru.f_2.D	    GmuStru.f_1.dc_1	GmuStru.f0.R        zeros(6)            zeros(6)            zeros(6)
            GmuStru.f_3.L	    GmuStru.f_2.dc1	    GmuStru.f_1.D	    GmuStru.f0.dc_1	    GmuStru.f1.R        zeros(6)            zeros(6)
            zeros(6)            GmuStru.f_2.L	    GmuStru.f_1.dc1	    GmuStru.f0.D        GmuStru.f1.dc_1	    GmuStru.f2.R        zeros(6)
            zeros(6)            zeros(6)            GmuStru.f_1.L	    GmuStru.f0.dc1	    GmuStru.f1.D        GmuStru.f2.dc_1	    GmuStru.f3.R
            zeros(6)            zeros(6)            zeros(6)            GmuStru.f0.L        GmuStru.f1.dc1	    GmuStru.f2.D        GmuStru.f3.dc_1
            zeros(6)            zeros(6)            zeros(6)            zeros(6)            GmuStru.f1.L        GmuStru.f2.dc1	    GmuStru.f3.D ];
    end

    A_hss = A_Tpz + C_Tpz * Gmxh*1 - F_shift; %% including CCSC control and possibly also AC current control
    B_hss = B_Tpz + C_Tpz * Gmuh*1; %% including controls with AC/DC voltage inputs

    H_hss = inv( s*eye(84) - A_hss ) * B_hss;
    %% By matrix inversing always consider the harmonic coupling
    if 1 %MMC_Y_PORT_NUM == 3
        %% Setting zero sequence ac current as 0 or eliminating it in the HSS matrix does not change the following result 
        H_hss_pnz = Tabc2pnz_hss_x * H_hss * Tpnz2abc_hss_u;
        %% Only extracting the ip,in,idc <-> vp,vn,vdc as follows can correctly integrate the CCSC impedance shaping effect (checked on 23.06.2023)
        Ymmc_hss_new =  H_hss_pnz([1:2 6     13:14 18   25:26 30   37:38 42   49:50 54   61:62 66   73:74 78], [1 2 6  7 8 12  13 14 18    19 20 24    25 26 30    31 32 36    37 38 42]);
        %Ymmc_hss_new =  H_hss_pnz([1:2 5     12:13 16   23:24 27   34:35 38   45:46 49   56:57 60   67:68 71], [1 2 6  7 8 12  13 14 18    19 20 24    25 26 30    31 32 36    37 38 42]);
        Zmmc_hss_pnz = inv(Ymmc_hss_new);
        Zmmc_pn_pu = [  -Zmmc_hss_pnz(4*3+1,4*3+1)   -Zmmc_hss_pnz(4*3+1,2*3+2)   Zmmc_hss_pnz(4*3+1,3*3+3)/3
                        -Zmmc_hss_pnz(2*3+2,4*3+1)   -Zmmc_hss_pnz(2*3+2,2*3+2)   Zmmc_hss_pnz(2*3+2,3*3+3)/3
                        -Zmmc_hss_pnz(3*3+3,4*3+1)	 -Zmmc_hss_pnz(3*3+3,2*3+2)	  Zmmc_hss_pnz(3*3+3,3*3+3)/3   ];% * Zbase;
        %% the minus sign assigned to H_hss_pnz can also be moved into Zmmc_pn_pu
    end
    Ymmc_pn_pu = inv(Zmmc_pn_pu); % in pu

    % if mmcSetup.CtrlMode == GFM_VF %|| VSC_CTRL_MODE == GFM_APC
    %     %% including terminal transformer into MMC impedance
    %     Ztr_pu = [Rxfo+Lxfo*(s+1j*wb)/wb   0; 0     Rxfo+Lxfo*(s-1j*wb)/wb ]; % * Zbase;
    %     Zmmc_pn_pu(1:2,1:2) = (Zmmc_pn_pu(1:2,1:2) + Ztr_pu); % * IS_Z_WITH_TRAFO);
    %     %Zmmc_pn = Zmmc_pn * (Vnom_prim/Vnom_sec)^2; %% 
    % end

    % 
    %     %% Simplified 2port AC modelling without considering DC dynamic
    % if 0 && MMC_Y_PORT_NUM == 2 %% for AC admittance modelling
    %     %% extracting frequency-shifted sequence impedances under modified sequence domain similar to rotor rotating dq frame
    %     %%          Zpp or Z11           Zpn or Z12          Znp or Z21               Znn or Z22
    %     Ymmc_hss = - H_hss([1:3     13:15   25:27   37:39   49:51   61:63   73:75], [1:3     7:9     13:15   19:21   25:27   31:33   37:39]);
    % 
    %     Zmmc_hss_pnz = Tabc2pnz_hss * inv(Ymmc_hss) * Tpnz2abc_hss; %% in real sequence domain
    %     Zmmc_pn_pu = [  Zmmc_hss_pnz(4*3+1,4*3+1)  Zmmc_hss_pnz(4*3+1,2*3+2)
    %                     Zmmc_hss_pnz(2*3+2,4*3+1)  Zmmc_hss_pnz(2*3+2,2*3+2)];% * Zbase; %% w1 rotating frame;
    %     Ymmc_pn_pu = inv(Zmmc_pn_pu);
    % end

    Zmmc_pu_buf(Cnt, :) = reshape(Zmmc_pn_pu.',1,[]);
    Ymmc_pu_buf(Cnt, :) = reshape(Ymmc_pn_pu.',1,[]);
    
    % %% Simplified case: Neglecting internal multi-frequency-coupling dynamics
    % for tt=1:1
    % if 0    
    %     Zrl_pn = [Rff+Lff*(s+1j*wb)/wb      0;      0       Rff+Lff*(s-1j*wb)/wb];
    %     Zrl_dq = [Rff+Lff*s/wb      -Lff;   Lff       Rff+Lff*s/wb];
    %     Gd = 1/(1+1.5*s*Ts);
    %     Hpll = CCCTRL.PLL.kp + CCCTRL.PLL.ki/s; Gpll = Hpll/(s+Hpll*U1d);
    %     Gpllu=[1  Gpll*U1q; 0  1-Gpll*U1d]; Gplli=[0  Gpll*I1q; 0  -Gpll*I1d]; Gpllm=[0  -Gpll*M1q; 0  Gpll*M1d];
    %     Hcc = Kp_Ireg + Ki_Ireg/s; Gcc = [ Hcc       Lff;  -Lff      Hcc];
    %     H_Pctrl = Kp_Preg + Ki_Preg/s; %% for positive sequence
    %     %% Q control with -1 coefficient for the PI control output
    %     H_Qctrl = -Kp_Qreg - Ki_Qreg/s; %% for positive sequence
    %     K2 = [0   1;    -1   0];
    %     G_PQu = [-H_Pctrl*i0_ac.'; -H_Qctrl*i0_ac.'*K2];
    %     G_PQi = [-H_Pctrl*u0_ac.'; H_Qctrl*u0_ac.'*K2];
    %     Zmmc_pn = (Zrl_pn + 0*Kpwm*Gcc/Udc0) * Zbase;
    %     H_LPF = (2*pi*Fn_filter)^2/(s^2+2*Zeta_filter*(2*pi*Fn_filter)*s+(2*pi*Fn_filter)^2);
    %     Zmmc_dq = inv( eye(2) - Gd*H_LPF*1.0*eye(2)) * ( Zrl_dq + Gcc*H_LPF ) * Zbase;
    %     Km=1; Vdc0=1;
    %     Zmmc_dq = inv( eye(2) - Km/Vdc0*Gd * ( 2*Gpllm + (Hcc*G_PQu+eye(2))*Gpllu*H_LPF+(Hcc*G_PQi-Gcc)*Gplli*H_LPF ) ) * ( Zrl_dq + Km/Vdc0*Gd * (Gcc-Hcc*G_PQi)*H_LPF ) * Zbase;
    %     d_theta0 = 0 - theta0;
    %     Tdq_l2g = [cos(d_theta0) sin(d_theta0); -sin(d_theta0) cos(d_theta0)];
    %     Zmmc_dq = Tdq_l2g * Zmmc_dq * inv(Tdq_l2g);
    %     Zmmc_pn = inv(Tpn2dq) * Zmmc_dq * Tpn2dq;         
    % end
    % end
    %% END
end
disp('Analytical Modelling Completed !')

YmmcStru.f_buf = f_buf;
YmmcStru.Zmmc_pu_buf = Zmmc_pu_buf;
YmmcStru.Ymmc_pu_buf = Ymmc_pu_buf;

%% if an impedance file does exist, save the impedance data for analysis
if ~isfile(IMPEDANCE_FILE)
    save(IMPEDANCE_FILE,"YmmcStru");
end
