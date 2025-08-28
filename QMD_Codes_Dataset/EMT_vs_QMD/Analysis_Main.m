%% EMT Simu (get Measurents for SS values Z responses) - Analytical Modelling - Validation (Using EMT Measurents)
addpath([pwd '/lib'])

%VSM = initVSM;
my_MMC_param

global SEQUENCE_NO ANALYTICAL_FLAG MEAS_FLAG
SEQUENCE_NO = 1; ANALYTICAL_FLAG = 0; MEAS_FLAG = 1;

%% General MMC settings
MMC.CtrlMode = GFL_APC; %% selection of control mode: GFM_VF = 0; GFL_APC = 1; GFL_DVC = 2; GFM_APC = 3; (values programmed in EMT models)
MMC.Ctrl.Modulation = 1; % 1 UdcN OR Udc0, 2 udc (or averaged udc??), 3, arm capacitor voltage sum
MMC.Ctrl.CCSC.Flag = 1; %% 0 for disabled, 1 for closed-loop control, 2 for open-loop control - not used anymore since 2024
MMC.Ctrl.ZSCC.Flag = 1; %% 0 for disabled, 1 for energy control, 2 for Idc active damping (active if CCSC_FLAG = 1)
MMC.Ctrl.IsPacFFinWt = 1; %% 1 denotes having power feedforward in energy control, 0 for not
%% Settings for GFL
MMC.Ctrl.IsCrossGFL = 0; %% 1 for cross control, 0 for classic control
MMC.Ctrl.IsVacGFL = 0; %% 1 for PV GFL mode, 0 for PQ GFL mode
%% Settings for GFM
MMC.Ctrl.IsCCinGFM = 0; %% 0 for without inner loop, 1 for with inner loop
MMC.Ctrl.IsVIinGFM = 1; %% 1 for having virtual impedance (VI), and 0 for not
MMC.Ctrl.IsPfDroopGFM = 0; %% droop is only enabled e.g. when f is not within [49.8 50.2] Hz, as implemented in EMT model
%% Setpoints
MMC.Set.Pref = 0.6; %1001e6/VSC1.pu.Sb; positive for feeding into AC-grid
MMC.Set.Qref = 0.0; %0*61e6/VSC1.pu.Sb;

%% create MMC ojects
%% !!! my_MMC_Impedance() only support MMC as input now
MMC1 = MMC;
mStru = my_MMC_Modeling(MMC);
%
MMC.CtrlMode = GFL_DVC;
MMC.Set.Pref = -0.6; %% power feeding into AC-grid as positive
%MMC.Ctrl.CCSC.Flag = 0;
MMC2 = MMC;
mStru2 = my_MMC_Modeling(MMC);

%% Component level analysis
INTERPOLATION_FLAG = 0;

% close %% will close all existing plots
if 1
    PLOT_Y_OR_Z_FLAG = 2; %% 1 for Y, 2 for Z
    mStru = mStru2;
    figure,
    %% Plot MMC analytical impedances/admittances
    plot_analytical_fun(mStru.f_buf,mStru.Zmmc_pu_buf,mStru.Ymmc_pu_buf,PLOT_Y_OR_Z_FLAG);
    %% Get and Plot measured impedances/admittances
    if INTERPOLATION_FLAG == 1
        f_buf = mStru.mMeas.Zmeas_pu_buf(:,1);
        f_buf_query = 0.1:0.1:100;
        append_part = 101:1:3000;
        f_buf_query(end+1:end+length(append_part)) = append_part;
        f_Zmeas_buf(:,1) = f_buf_query.';
        f_Ymeas_buf(:,1) = f_buf_query.';
        for curNo=2:1:10
            f_Zmeas_buf(:,curNo) = interp1(f_buf,mStru.mMeas.Zmeas_pu_buf(:,curNo),f_buf_query);
            f_Ymeas_buf(:,curNo) = interp1(f_buf,mStru.mMeas.Ymeas_pu_buf(:,curNo),f_buf_query);
        end
        plot_meas_fun(f_Zmeas_buf,f_Ymeas_buf,PLOT_Y_OR_Z_FLAG);
    else
        plot_meas_fun(mStru.mMeas.Zmeas_pu_buf,mStru.mMeas.Ymeas_pu_buf,PLOT_Y_OR_Z_FLAG);
    end
    return,
end
%my_RGA_fun(mStru.f_buf, mStru.Ykk_pu_buf);

if ANALYTICAL_FLAG == 0 && MEAS_FLAG == 1
    f_buf = mStru.mMeas.Ymeas_pu_buf(:, 1);
    Ymmc1_pu_buf = mStru.mMeas.Ymeas_pu_buf(:, 2:end);
    Ymmc2_pu_buf = mStru2.mMeas.Ymeas_pu_buf(:, 2:end);
    if INTERPOLATION_FLAG == 1
    %% interpolation, added in 06.2025
        if Hmax <= 2
            f_buf_query = Hmin*50:0.1:Hmax*50;
        elseif Hmin < 2
            f_buf_query = Hmin:0.1:100;
            append_part = 101:1:Hmax*50;
            f_buf_query(end+1:end+length(append_part)) = append_part;
        else % if Hmin >= 2
            f_buf_query = Hmin*50:1:Hmax*50;
        end
        for curNo=1:1:9
            Ymmc1_finer(:,curNo) = interp1(f_buf,Ymmc1_pu_buf(:,curNo),f_buf_query);
            Ymmc2_finer(:,curNo) = interp1(f_buf,Ymmc2_pu_buf(:,curNo),f_buf_query);
        end
        f_buf = f_buf_query;
        Ymmc1_pu_buf = Ymmc1_finer;
        Ymmc2_pu_buf = Ymmc2_finer;
    end
    %% end of interpolation
else
    f_buf = mStru.f_buf;
    Ymmc1_pu_buf = mStru.Ymmc_pu_buf;
    Ymmc2_pu_buf = mStru2.Ymmc_pu_buf;
end

%Pac1_vec = [3 2.6 2.4 2.3 2.2]*Pnom;
Pac1_vec = [2.5]*Pnom;
Pac2_vec = [6]*Pnom;
for cur_no=1:1:length(Pac1_vec)
%% -------------------  AC/DC grid (excluding MMCs) and system modeling  -------------------------
%% Grids AC1 and AC2
% Pac1 = Psc;
% Pac2 = Psc2;
%Pac1 = Pnom * 2.5;
Pac1 = Pac1_vec(cur_no);
Pac2 = Pnom * 6;
%Pac2 = Pac2_vec(cur_no);
%% Simulink calculates the L-based SCL, resulting Xgrid = Zg_mag
Zg_mag = Vnom_prim^2/(Pac1*1); Xgrid = Zg_mag;% * X_R/sqrt(X_R^2+1^2); %% abs
%Xgrid = Zg_mag * X_R/sqrt(X_R^2+1^2); %% abs
Rgrid = Xgrid/X_R; Lgrid = Xgrid/wb; %% abs
Rg_pu = 1*Rgrid*(Vnom_sec/Vnom_prim)^2/Zbase;
Lg_pu = 1*Lgrid*(Vnom_sec/Vnom_prim)^2/(Zbase/wb); %% pu under system base
Rg2_pu = Rg_pu*(Pac1/Pac2);
Lg2_pu = Lg_pu*(Pac1/Pac2); %% pu under system base

for k=1:1:length(f_buf)
    
    f=f_buf(k); % Hz
    w=2*pi*f;
    s=j*w;
    %% AC grid
    Zgac_pu = [Rxfo+Rg_pu+(s+1j*wb)*(Lxfo+Lg_pu)/wb            0
                0                   Rxfo+Rg_pu+(s-1j*wb)*(Lxfo+Lg_pu)/wb];
    Zgac2_pu = [Rxfo+Rg2_pu+(s+1j*wb)*(Lxfo+Lg2_pu)/wb            0
                0                   Rxfo+Rg2_pu+(s-1j*wb)*(Lxfo+Lg2_pu)/wb];
    Ygac_pu = inv(Zgac_pu);
    Ygac2_pu = inv(Zgac2_pu);

    %% DC grid: sinlge-pole mapped to double-pole through dividing Ydc by 2 or equivalently as follows (induced by replacing vdc/2 with vdc)
    if 1 %% used for MMC impedance validation
        Zgdc_pu = 2*(R_cable+s*L_cable)/Zbase; %% + 2*1/(s*C_PM/Nb_PM*3 + 1/(s*Ldc_r))/Zbase;
        Rload = (640)^2/600; %% Ohm, as the sum of Rloads +/- = (320/2)^2/300; %% Ohm
        if MMC1.CtrlMode == GFL_DVC
            Zgdc_pu = Rload/Zbase + Zgdc_pu;
        end
        Ygdc_pu = 1/Zgdc_pu; %% Zbase is used, which should be replaced by Zdcb in later studies
    end
    
    YpDC = s*Cdc_perKm*Length_DC1/2/Ybase / 2;
    ZsDC = (Rdc_perKm+s*Ldc_perKm)*Length_DC1/Zbase * 2;
    Zdc_r = s*Ldc_r/Zbase * 2;

    [~, kxk] = size(Ymmc1_pu_buf);
    k_dim = sqrt(kxk);
    Ymmc1_pu = reshape(Ymmc1_pu_buf(k, :),[k_dim,k_dim]).';
    Ymmc2_pu = reshape(Ymmc2_pu_buf(k, :),[k_dim,k_dim]).';
    

    %% Ygac1 and Ygac2 should be modelled
    %%                   AC1                         DC1                        DC2                             AC2
    Ykk_pu = [Ymmc1_pu(1:2,1:2)+Ygac_pu      Ymmc1_pu(1:2,3)                    zeros(2,1)                      zeros(2)
              Ymmc1_pu(3,1:2)       1/(1/Ymmc1_pu(3,3)+Zdc_r)+YpDC+1/ZsDC       -1/ZsDC                         zeros(1,2)
              zeros(1,2)                        -1/ZsDC             1/(1/Ymmc2_pu(3,3)+Zdc_r)+YpDC+1/ZsDC       Ymmc2_pu(3,1:2)
              zeros(2)                          zeros(2,1)                       Ymmc2_pu(1:2,3)                Ymmc2_pu(1:2,1:2)+Ygac2_pu];

    %%                   AC1                         DC1                        DC2                             AC2
    Ykk_pu0 = [Ymmc1_pu(1:2,1:2)+Ygac_pu      Ymmc1_pu(1:2,3)                    zeros(2,1)
              Ymmc1_pu(3,1:2)       1/(1/Ymmc1_pu(3,3)+Zdc_r)+YpDC+1/ZsDC       -1/ZsDC 
              zeros(1,2)                        -1/ZsDC             1/(0/Ymmc2_pu(3,3)+Zdc_r)+YpDC+1/ZsDC];
    %% When simplifying MMC2 as an ideal DC voltage source
    Ygdc_pu = 1 / ( Zdc_r+1/(YpDC+1/(ZsDC+0/(YpDC+1/Zdc_r))) );
    %Ygdc_pu = 1 / ZsDC;

    %% Save frequency-sweep results in buffers for further analysis
    Ygac_pu_buf(k, :)= reshape(Ygac_pu.',1,[]);
    Ygac2_pu_buf(k, :)= reshape(Ygac2_pu.',1,[]);
    Ygdc_pu_buf(k) = Ygdc_pu;
    Ykk_pu_buf(k, :)= reshape(Ykk_pu.',1,[]);
end

%% System Modelling
for k=1:1:length(f_buf)
    f=f_buf(k); % Hz
    w=2*pi*f;
    s=j*w;

    [~, kxk] = size(Ymmc1_pu_buf);
    k_dim = sqrt(kxk);
    Ymmc1_pndc = reshape(Ymmc1_pu_buf(k, :),[k_dim,k_dim]).';
    Ymmc2_pndc = reshape(Ymmc2_pu_buf(k, :),[k_dim,k_dim]).';

    % if MMC_Y_PORT_NUM == 2
    %     %% AC equivalent considering DC dynamic
    %     y_dc2ac = Ymmc_pn_pu(1:2,3);
    %     y_ac2dc = Ymmc_pn_pu(3,1:2);
    %     Ymmc_pn_pu = Ymmc_pn_pu(1:2,1:2) - y_dc2ac * inv(Ygdc_pu+Ymmc_pn_pu(3,3)) * y_ac2dc;
    %     Zmmc_pn_pu = inv(Ymmc_pn_pu);
    % end
    % 
    % if MMC_Y_PORT_NUM == 2
    %     %Ymmc_pn(1,2) = 0; Ymmc_pn(2,1) = 0;
    %     Ykk_pu = Ymmc_pn_pu + 1*Ygac_pu;
    % elseif MMC_Y_PORT_NUM == 3
    %     Ykk_pu = [   Ymmc_pn_pu(1:2,1:2)+Ygac_pu*1      Ymmc_pn_pu(1:2,3)
    %                  Ymmc_pn_pu(3,1:2)                  Ymmc_pn_pu(3,3)+Ygdc_pu*1 ];
    % end
end

%% the algorithm used for SISO impedances still need to be further verified - 05-2024
if 0
my_SISO_fun(f_buf,mStru.Zmmc_pu_buf,Ygac_pu_buf,Ygdc_pu_buf,mStru.mMeas);
%my_SISO_fun(f_buf,mStru2.Zmmc_pu_buf,Ygac2_pu_buf,Ygdc_pu_buf,mStru2.mMeas);
return,
end
%% --------------------------------------------------------------------
RMA_CURVE_FLAG = 0;
RMAout = my_RMA_fun(f_buf,Ykk_pu_buf,RMA_CURVE_FLAG);
if RMA_CURVE_FLAG == 1
    continue,
end
fk_RMA = RMAout(:,5);
D_RMA = RMAout(:,6);
%fk_RMA2pole = RMAout(:,7);
semilogx(fk_RMA, D_RMA,'rx'), hold on, grid on
xlabel('{\it f} in Hz'), ylabel(['\zeta in pu']);
title('P2P-HVDC test system')
legend('Modes (Z-Curves)')
%my_ACDC_int(f_buf, Zmmc_pu_buf, Ygac_pu_buf, Ygdc_pu_buf, mMeas);
%det_stability_fun(f_buf,Ykk_sys_buf);
%% --------------------------------------------------------------------
end