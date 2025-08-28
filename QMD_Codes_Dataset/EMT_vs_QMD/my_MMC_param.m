%% Definition of control modes
GFM_VF = 0; GFL_APC = 1; GFL_DVC = 2; GFM_APC = 3;
%
Fnom=50; Vnom_prim=400e3;%400e3;       % Nominal primary voltage (V)
Pnom=1000e6; Vnom_sec=320e3; Vnom_dc=640e3;
%% GLOBAL TEST SETTINGS
Ts_Power=20e-6 * 5;
Ts_Control=Ts_Power*2; Ts=Ts_Control;
PhaseDelay = 1.0*Ts_Control*2*pi*50;%% Set a large delay in EMT simulation can induce instability (e.g. in case of large Ts_Control)
% Analytical modeling setup
Hmin = 1/50; Hstep = Hmin/100; Hmax = 0.6; % pu

IS_VNO_INCLUDED = 0;
MMC_Y_PORT_NUM = 3; %% 3 for hybrid AC-DC admittance matrix, 2 for AC admittance matrix, 1 for DC admittance
%GFL_Iref_IS_PQref = 1; %% setting types of PQ controllers when implemented
SIMU_Vinj_FLAG = 1;

fb=50; wb=2*pi*fb;%base angular frequency
%% --------------------------------------
%% setting grid conditions
f0=50/fb;
U_RMS_pu=1; %% furtherly reducing it for current-limiting stability analysis
%% end here
w0=f0;%wb*1.0000;
%% grid frequency that may deviate from 50 Hz Used in the Injection Toolbox !!!
f1_Hz_Vinj=f0*fb; 
%% initVSM
VSM.Ta  = 10*1;   % Ta    Time required to spin up in [s]
VSM.kpf = 20*1;   % kpf = deltaP/deltaf  in  [pu/pu]       (reasonable: 2.5)
%% dividing kdf by 10 greatly reduce damping, which can be used to create instability
VSM.kdf = 0.125*1;  % kdf = deltaPhi/deltaf in [rad/(rad/s)]
VSM.f0  = 50;   % f0  = nominal frequency  [Hz]
VSM.ku  = 0.1*1;  % ku  = deltaU/deltaQ  in  [pu/pu]
VSM.Tu  = 0.3;    % time constant of Q/U-loop (AVR equivalent Time Constant)
VSM.VC.tau_rv = 1/(20*2*pi);
VSM.VC.rv = 0.13*1;
VSM.VC.xv = 0.11*1;
%VSM.PLL = CCCTRL.PLL;
VSM.VTFF.OmLead = 2*pi*50*0;
VSM.VTFF.OmLag  = 2*pi*1100;
VSM.IPF.Deadband = 0.2/50;
%% negative by Q-control with negative gain, positive by by Q-control with positive gain
VSM.k_qu  = -0.2*1; %% its variation influences f<1Hz dynamics
VSM.kQ = -1*1; %% its variation influences f<1Hz dynamics
%% cut-off frequency for PQ_LPF, which is actually not needed as the control structure of VSM from Daniel und Selfsync from Unruh are not LPF-based
VSM.w_c = 2*pi*5e3; %% thus here setting a very large value, 5e3Hz,as default, indicating disabled
%
%% ---------------- AC Grid Setting ----------------------
Psc = Pnom * 6;  % Short circuit power (MVA), used for measuring impedance, in Simulink only representing X-part
X_R = 10; % X/R ratio, in EHV grid generally > 10: as for Tr X_R >= 50, for lines X_R >= 10 (headline/cable ) 
%Psc2 = Pnom * 20;
%% -----------------DC Cable data ----------------------
R_cable = 0.5 * 5;    % ohm
L_cable = 15e-3 * 20;   % (H)
%% -------------------------------------------------------
% Sequencer timing:
Tbrk1_On=0.1;                 % Closing time of breaker 1 (converter energizing)
Tbrk2_On=1.0;                 % Closing time (s) of breaker 2 (across start-up resistor)
%
Tdeblock=1.5;                 % Converter de-block time (s)
Ton_VDCreg=1.5;               % VDC regulator turn-on time (s) - VDC Regulation
Tramping_Vdc_ref=2;           % Start time Vdc_ref ramping to nominal (s)
Slope_Vdc_ref=Vnom_dc/2;      % Slope Vdc_ref ramping (V/s)
%
Ton_PQreg=4;                  % Preg & Qreg regulators turn-on time (s), Inherited also by Energy control and other control modes
Tramping_Pref=Ton_PQreg+0.2;  % Start time Pref ramping(s)
Slope_Pref=0.5;               % Slope Pref ramping (V/s)
Tramping_Qref=Ton_PQreg+2;  % Start time Qref ramping(s)
Slope_Qref=0.5;               % Slope Pref ramping (V/s)
Ton_CROSSreg=Ton_PQreg+3;
%
Ton_Converter2=4;             % Converter 2 equivalent switched-on time (s)
%
Nb_PM = 400;                      % Number of power module per arm
C_PM = 10e-3; % Power module capacitor (F) 
% Energy in kJ/MVA
W_kJ_MVA= 0.5 * C_PM * (Vnom_dc/Nb_PM)^2 * Nb_PM * 6 / (Pnom/1e6)/1e3;
Wt_ref_pu= 0.5 * C_PM * (Vnom_dc/Nb_PM)^2 * Nb_PM * 6 / Pnom; %% Wt in pu is also defined as electrostatic constant Hdc (typically varying from 5ms to 40ms)
Vc0_PM=0;                     % Capacitors initial voltage (V)
% Transformer impedance
Lxfo= 0.16;         % Total Leakage inductance (pu), default setting is 0.12 (it is 0.18 in INELFE link)
Rxfo= 0.002;        % default value is 0.003, Total winding resistance (pu)
Zbase= Vnom_sec^2/Pnom; Ybase=1/Zbase;
%% default setting for Larm_pu is 0.15 pu (same as in INELFE link)
Larm_pu=0.15; Larm=Larm_pu*(Zbase/wb);
Rarm_pu=Larm_pu/100; Rarm=Rarm_pu*Zbase; %% Lpu=L/(Zb/wb)=L*wb/Zb=Xpu
%% DC base to be checked
Sbase=Pnom;
Ubase=sqrt(2/3)*Vnom_sec; 
Ibase=2*Sbase/(3*Ubase); %Ib = Sb / Ub;
Kpwm = 0.5;
%% As all control variables are pu values obtained using AC base system, DC base should be carefully used.
%%  Udc0=Vnom_dc/Udcb;%in pu value
%% DC BASE not used anymore. The whole system then employs the AC base
Udcb=2*Ubase; Idcb=Sbase/Udcb;%=3/4*Ib; %%% Sb=1.5*Ub*Ib=Udcb*Idcb;
Udc0=Vnom_dc/Udcb;%in pu value
Zdcb=Udcb/Idcb; %%% 8/3*Zb;
Ydcb=1/Zdcb; Cdcb = 1/Zdcb/wb; Ldcb = Zdcb/wb;
%% Multiple SM capacitors are not like the single C in TL-VSC DC side, instead every arm has an equivalent Carm (thus in abc phases) 
Cbase = 1/Zbase/wb; %% defined in the way of Yc_pu = s*Cpu/wb like Zl_pu = s*Lpu/wb
Csm_pu = C_PM / Cbase; %% in per unit value
%
%% dq and Vdc measurement filter cut-off frequency:
Fn_filter=1000; % in Hz
Zeta_filter=1;

%% PLL: slow PLL can induce startup instability
CCCTRL.PLL.T  = 0.01; % default 0.04 s
CCCTRL.PLL.v0 = 1;
Kp_PLL = 2/CCCTRL.PLL.v0/CCCTRL.PLL.T;
Ki_PLL = 1/(CCCTRL.PLL.v0*CCCTRL.PLL.T^2);
%
%% Active power regulator (Preg) : tau=100 ms
Kp_Preg= 0.5/3*1;                % Proportional gain
Ki_Preg= 1.0*10;                  % Integral gain
Limits_Preg = [ 2, -2 ] ;   % Output (Vdc_ref) Upper/Lower limits (pu)
%
%% Reactive power regulator (Qreg): tau=100ms
Kp_Qreg= 0.5/3*1;                % Proportional gain
Ki_Qreg= 1.0*10;                  % Integral gain
Limits_Qreg = 2*[ 0.42, -0.42 ]; % Output (Iq_ref) Upper/Lower limit (pu)
%
%% VDC regulator (VDCreg): tau=100ms, high bandwidth can make impedance extraction from EMT simulation difficult
Limits_VDCreg= [ 2.0  -2.0];   % Output Idref [Upper Lower] limits (pu)
Kp_VDCreg=1;                   % Proportional gain
Ki_VDCreg=10;                 % Integral gain
%% the following larger parameters have better performance (02.08.2024)
% Kp_VDCreg=4;                   % Proportional gain
% Ki_VDCreg=100;

%% VAC regulator for VF mode
Kp_VACreg = 0.01;
Ki_VACreg = 10;
Limits_VACreg= [ 1.3  -1.3];
%
%% VAC for GFL PV mode: tau=100ms 
Kp_VAC_GFL = 0.5;
Ki_VAC_GFL = 50;
%
%% AC Current regulator (Ireg):tau=1ms
Limits_Ireg= [ 2.0  -2.0];     % Output Vdq_conv [Upper Lower] limits (pu)
CCCTRL.CC.design.SCR = 4;
CCCTRL.CC.design.XR = 10;
CCCTRL.CC.design.TCC = 2e-3; % default 0.001 s
CCCTRL.CC.design.LCL.PU.L2 = 1/CCCTRL.CC.design.SCR;% pu 
CCCTRL.CC.design.LCL.PU.R2 = 1/CCCTRL.CC.design.XR/CCCTRL.CC.design.SCR;% pu 
CCCTRL.CC.kp = (Larm_pu/2 + CCCTRL.CC.design.LCL.PU.L2) / wb / CCCTRL.CC.design.TCC  ;
CCCTRL.CC.ki = (Rarm_pu/2 + CCCTRL.CC.design.LCL.PU.R2) / CCCTRL.CC.design.TCC ;
Kp_Ireg = CCCTRL.CC.kp; %% 
Ki_Ireg = CCCTRL.CC.ki; %% 
%
%% ZSCC controller (icir_z):tau=1ms, updated on 17.09.2024 (Influencing LF passivity)
Kp_Izscc = 10;
Ki_Izscc = 100;

%% CCSC (icir_2dq):tau=5ms, Updated on 17.09.2024 -tuned by checking the time of arriving 0.66*Step_value
Kp_Iccsc = 0.21;
Ki_Iccsc = 5.15;

%% Energy Controller:tau=50ms, tuning criterion is 5*tau to reach 98% steady-state value for a 2nd-order system response
Kp_WtReg = 40;
Ki_WtReg = Kp_WtReg / 0.025;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ZSCC active damping
Rad_ZSCC = 0.5;
Tad_ZSCC = 1/(2*pi*5);

%% Feedforward coefficients:
Lff=Larm_pu/2;
Rff= Rarm_pu/2;
%
%% Added GFM inner loop V Control
Kpv_GFM = 1;
%
%% GFL QU droop control
CCCTRL.QU.Deadband = 0.03; %% nach VDE AR-N 4105 for Q(U) control
CCCTRL.QU.k = -10*1; %with grid-following (GFL) control ; %% Grid Code schauen und überprüfen, ob es positiv oder negativ sein soll
CCCTRL.QU.T = 6/3; %% nach VDE AR-N 4105 zu parametrieren
%
%% Start up settings
P_Ld1 = Psc/30;   % load (primary bus) (MW)
R_startup = 400;   % Startup resistance (Ohm)
%
%
%% DC smoothing reactor
Ldc_r=60e-3; % H
Length_DC1=200; % km
%% selected parameters for DC line/cable +/-320kV
Rdc_perKm=1.1e-3; % Ohms
Ldc_perKm=2e-4; % H
Cdc_perKm=1.5e-8; % F
if 0 %% DC line, parameter source: ???
Rdc_perKm=1.39e-2; % Ohms
Ldc_perKm=1.59e-4; % H
Cdc_perKm=1.31e-7; % F
end
if 1 %% paramter set from CIGRE TB604 for DC cable +/-400kV
    Rdc_perKm=1.1e-2; % Ohms
    Ldc_perKm=2.6e-3; % H
    Cdc_perKm=1.91e-8; % F
end
if 0 %% Parameters from the Book "Analysis and Mitigation of Broadband Oscillation in Renewable Energy Generation and AC/DC Transmission Systems", W.Wang,2023.
Rdc_perKm=3.3e-3; % Ohms
Ldc_perKm=2e-4; % H
Cdc_perKm=1.5e-8; % F
end
%% AC Cables: parameters from the EHV-Grid
Length_AC1 = Length_DC1; % km
Rac_perKm=2.5e-2; % Ohms
Lac_perKm=7.96e-4; % H, practicability to be checked
Cdc_perKm=1.4e-8; % F
%
%% Including Transformation matrices
Tabc2pn = [1    exp(j*2*pi/3)   exp(j*4*pi/3);  1    exp(j*4*pi/3)  exp(j*2*pi/3)]/3;
Tpn2abc = [1    1; exp(j*4*pi/3)   exp(j*2*pi/3);  exp(j*2*pi/3)  exp(j*4*pi/3)];
Tp_abc = [1    exp(-j*2*pi/3)   exp(j*2*pi/3)].';
Tn_abc = [1    exp(j*2*pi/3)   exp(-j*2*pi/3)].';

Tabc2pnz = [1    exp(j*2*pi/3)   exp(j*4*pi/3);  1    exp(j*4*pi/3)  exp(j*2*pi/3); 1   1   1]/3;
Tabc2pn0 = [1    exp(j*2*pi/3)   exp(j*4*pi/3);  1    exp(j*4*pi/3)  exp(j*2*pi/3); 0   0   0]/3;
Tpnz2abc = [1    1  1; exp(j*4*pi/3)   exp(j*2*pi/3)    1;  exp(j*2*pi/3)  exp(j*4*pi/3)    1];

Tabc2z = [1   1   1]/3; Tz2abc = [1; 1; 1];

Tpn2dq = [1    1;  -j     j]; %% inversed matrix is Tdq2pn = [1    j;  1     -j]/2;
Tdq2pn = inv(Tpn2dq);