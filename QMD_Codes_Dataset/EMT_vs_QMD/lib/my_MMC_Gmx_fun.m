function GmxStru = my_MMC_Gmx_fun(IN,iMMC)

my_MMC_param

%% Transfer function from x_abc to m_abc

%% m_abc = [mac     mdc]
%% x_abc = [iac_abc   icir*ones(1,3)    uc_abc_cm   uc_abc_dm].' 
%% u_abc = [uac_abc     udc*ones(1,3)].'
wb = 2 * pi * 50;
w1 = wb;

H_offset=1e-6; %% pu of harmonic order

u0_ac = IN.u0_ac;
i0_ac = IN.i0_ac;
m0_ac = IN.m0_ac;
theta0 = IN.theta0;
m0_ccsc = IN.m0_cir;
i0_cir2 = IN.i0_cir2;
Uc_sum0_f2dq = IN.uc_sum0_f2dq;
Uc_sum0_f0z = IN.uc_sum0_f0z;
Uc_diff0_f1dq = IN.uc_diff0_f1dq;

U1d=u0_ac(1); U1q=u0_ac(2); %%% Us_vec=U1d+jU1q, measured values
I1d=i0_ac(1); I1q=i0_ac(2); %%% Ig1_vec=Ig1d+jIg1q, measured values
M1d=m0_ac(1); M1q=m0_ac(2);
U0=sqrt(U1d^2+U1q^2);
P0 = i0_ac.'*u0_ac;

h_vec = [-3 -2 -1 0 1 2 3];

for k=1:1:length(h_vec)
    h = h_vec(k);
    s = IN.s;
    if s+2j*w1+j*h*w1 == 0 || s-2j*w1+j*h*w1 == 0
        s = s-1j*H_offset*wb;
    end
    
    %% CCSC control from icir to mdc
    Hccsc = Kp_Iccsc + Ki_Iccsc/(s+j*h*w1);
    %% frequency-shifted controller transfer function
    Hccsc1 = Kp_Iccsc + Ki_Iccsc/(s+j*h*w1-(-2j*w1)); %% for positive sequence
    Hccsc2 = Kp_Iccsc + Ki_Iccsc/(s+j*h*w1+(-2j*w1)); %% por negative sequence
    
    % Gd = 1/(1+1.5*s*Ts); %% inaccurate by large delay
    % Gd = exp(-1.5*s*Ts);
    % Gd1 = exp(-1.5*(s+j*h*w1-j*w1)*Ts);
    % Gd2 = exp(-1.5*(s+j*h*w1+j*w1)*Ts);
    Gd = exp(-0.5*s*Ts)*(1-exp(-1.5*s*Ts))/(1.5*s*Ts);
    K2 = [0   1;    -1   0];
    %% set H_LPF=1 when no LPF is implemented for u,i,p,q measurements
    H_LPF = (2*pi*Fn_filter)^2/((s+j*h*w1)^2+2*Zeta_filter*(2*pi*Fn_filter)*(s+j*h*w1)+(2*pi*Fn_filter)^2);
    H_LPF1 = (2*pi*Fn_filter)^2/((s+j*h*w1-j*w1)^2+2*Zeta_filter*(2*pi*Fn_filter)*(s+j*h*w1-j*w1)+(2*pi*Fn_filter)^2);
    H_LPF2 = (2*pi*Fn_filter)^2/((s+j*h*w1+j*w1)^2+2*Zeta_filter*(2*pi*Fn_filter)*(s+j*h*w1+j*w1)+(2*pi*Fn_filter)^2);
    % C_LPF1 = (2*pi*Fn_filter)^2/((s+j*h*w1-(-2j*w1))^2+2*Zeta_filter*(2*pi*Fn_filter)*(s+j*h*w1-(-2j*w1))+(2*pi*Fn_filter)^2);
    % C_LPF2 = (2*pi*Fn_filter)^2/((s+j*h*w1+(-2j*w1))^2+2*Zeta_filter*(2*pi*Fn_filter)*(s+j*h*w1+(-2j*w1))+(2*pi*Fn_filter)^2);
    I_LPF1 = 1; I_LPF2 = 1;
    C_LPF1 = 1; C_LPF2 = 1; C_LPF = 1; 

    Q_LPF1 = (2*pi*Fn_filter/10)^2/((s+j*h*w1-j*w1)^2+2*Zeta_filter*(2*pi*Fn_filter/10)*(s+j*h*w1-j*w1)+(2*pi*Fn_filter/10)^2);
    Q_LPF2 = (2*pi*Fn_filter/10)^2/((s+j*h*w1+j*w1)^2+2*Zeta_filter*(2*pi*Fn_filter/10)*(s+j*h*w1+j*w1)+(2*pi*Fn_filter/10)^2);
    
    if iMMC.Ctrl.CCSC.Flag == 1
        %% impact from Gd_1 and Gd_2 proved negligible as compared to Gd 
        %Gd_1 = 1/(1+1.5*(s+j*h*w1-(-2j*w1))*Ts);
        %Gd_2 = 1/(1+1.5*(s+j*h*w1+(-2j*w1))*Ts);
        %% non-diagonal terms of the following matrix is obtained using the p/n relationships 
        Gccsc1 = [  Hccsc1      -2*Lff
                    +2*Lff      Hccsc1  ] * C_LPF1 * Gd / Udc0;
        %% The following transfer function is formulated based on the CCSC control structure        
        Gccsc2 = [  Hccsc2      -2*Lff
                    2*Lff      Hccsc2  ] * C_LPF2 * Gd / Udc0;
                
        %% performing the linear transformation from dq-frame to rotating modified sequence domain
        Gccsc_mpn1 = inv(Tpn2dq) * Gccsc1 * Tpn2dq;
        Gccsc_mpn2 = inv(Tpn2dq) * Gccsc2 * Tpn2dq;
        %% adding the static-angle rotation does not show any impedance shaping effect: because of no sequence coupling???
        
        theta0_f2 = -theta0*2;
        Tframe_f2 = [cos(0-theta0_f2) sin(0-theta0_f2); -sin(0-theta0_f2) cos(0-theta0_f2)];
        Tframe_f1 = [cos(0-theta0) sin(0-theta0); -sin(0-theta0) cos(0-theta0)];
        if 1
        %% Does Tdq2pn also apply to the rotation from -2dq rotating frame to sequence domain???
        Gccsc_mpn1 = inv(Tpn2dq) * Tframe_f2 * Gccsc1 * inv(Tframe_f2) * Tpn2dq;
        Gccsc_mpn2 = inv(Tpn2dq) * Tframe_f2 * Gccsc2 * inv(Tframe_f2) * Tpn2dq;
        end
        %% Frequency coupling effect from CCSC was found negligible
        %% frequency shifting already included in the frequency-shifted controller modelling
        Gccsc_pn = [Gccsc_mpn1(1,1)     1*Gccsc_mpn2(1,2)
                    1*Gccsc_mpn1(2,1)     Gccsc_mpn2(2,2)];
    
        Gccsc_abc = Tpn2abc * Gccsc_pn * Tabc2pn;
    
        %%% above modelling can be simplified to
        % if h == 2
        %     Gccsc_abc = 1 * (Kp_Iccsc + Ki_Iccsc/s - 2j*Lff) * C_LPF * Gd / Udc0 * eye(3);
        % elseif h == -2
        %     Gccsc_abc = 1 * (Kp_Iccsc + Ki_Iccsc/s + 2j*Lff) * C_LPF * Gd / Udc0 * eye(3);
        % else
        %     Gccsc_abc = zeros(3);
        % end
    
        if iMMC.Ctrl.ZSCC.Flag == 1
            %% add the single-loop ZSCC active damping control
            Hicir_z = Kp_Izscc + Ki_Izscc/(s+j*h*w1);
            K3 = diag([0 0 1]);
            Gwt_icir_z = Tpnz2abc * Hicir_z * C_LPF * K3 * Tabc2pnz;
            %% adding modulation compensation
            Gwt_icir_z = Gwt_icir_z / Udc0;
            Gccsc_abc = Gccsc_abc + Gwt_icir_z;
            %% Add the Wt response to Uc perturbations
            Hwt = Kp_WtReg + Ki_WtReg/(s+j*h*w1);
            %Gwt_uc = Hicir_z * Hwt * H_LPF * Csm_pu / wb / Nb_PM / 3 / Udc0_fUb;
            Gwt_uc = Hicir_z * Hwt * C_LPF * Csm_pu / wb / Nb_PM / 3 / Udc0 * Idcb / Ibase;
            %% adding modulation compensation
            Gwt_uc = Gwt_uc / Udc0;
            Gwt_uc_sum_z = Tpnz2abc * Gwt_uc * Uc_sum0_f0z * K3 * Tabc2pnz;
            %% effects of dq perturbation terms are negligible
            k3_vec = [0 0 1].';
            Gwt_uc_sum_dq = Tpnz2abc * Gwt_uc * k3_vec * Uc_sum0_f2dq.' * Tpn2dq * Tabc2pn;
            %Gwt_uc_sum = Gwt_uc_sum_dq + Gwt_uc_sum_z;
            Gwt_uc_diff = Tpnz2abc * Gwt_uc * k3_vec * Uc_diff0_f1dq.' * Tpn2dq * Tabc2pn;
        elseif iMMC.Ctrl.ZSCC.Flag == 2
            %% add the single-loop ZSCC active damping control
            Hzscc = Rad_ZSCC*(s+j*h*w1)/((s+j*h*w1)+1/Tad_ZSCC);
            K3 = diag([0 0 1]);
            Gmdc_icir_z = Tpnz2abc * C_LPF1 * Hzscc * K3 * Tabc2pnz / Udc0;
            Gccsc_abc = Gccsc_abc + Gmdc_icir_z;
            Gwt_uc_sum_dq = zeros(3);
            Gwt_uc_sum_z = zeros(3);
            Gwt_uc_diff = zeros(3);
        else %if ZSCC_FLAG == 0
            Gwt_uc_sum_dq = zeros(3);
            Gwt_uc_sum_z = zeros(3);
            Gwt_uc_diff = zeros(3);
        end
        
        %% choose if considering angle dynamics for CCSC, mainly from GFM synchronisation control
        if iMMC.CtrlMode ~= GFM_APC
            Gccsc_iac_3 = zeros(3);
            Gccsc_iac3 = zeros(3);
        else
            %% VSM angle dynamics (proved negligible)
            Hccsc1_f1 = Kp_Ireg + Ki_Ireg/(s+j*h*w1-(1j*w1)); %% for positive sequence
            Hccsc2_f1 = Kp_Ireg + Ki_Ireg/(s+j*h*w1+(1j*w1)); %% por negative sequence
        
            %% non-diagonal terms of the following matrix is obtained using the p/n relationships 
            Gccsc1_f1 = [  Hccsc1_f1      -2*Lff
                            +2*Lff      Hccsc1_f1  ] * C_LPF1 * Gd / Udc0;
            %% The following transfer function is formulated based on the CCSC control structure        
            Gccsc2_f1 = [  Hccsc2_f1      -2*Lff
                            2*Lff      Hccsc2_f1  ] * C_LPF2 * Gd / Udc0;
            % if abs(f0-1)<=0.2/50 
            %     VSM.k_pf=0;
            % else
            %     VSM.k_pf = VSM.kpf;
            % end
            LPF_INF = 1;
            Gsync1=-1*wb/(s+j*h*w1-j*w1)*(1+VSM.kdf*(s+j*h*w1-j*w1))/(VSM.Ta*(s+j*h*w1-j*w1)+iMMC.Ctrl.IsPfDroopGFM*VSM.kpf*iMMC.Set.Pref)*LPF_INF;
            Gsync2=-1*wb/(s+j*h*w1+j*w1)*(1+VSM.kdf*(s+j*h*w1+j*w1))/(VSM.Ta*(s+j*h*w1+j*w1)+iMMC.Ctrl.IsPfDroopGFM*VSM.kpf*iMMC.Set.Pref)*LPF_INF;
            %% when the active power is calculated using dq variables in VSM dq-frame
            G_SYNC1=Gsync1/(1-Gsync1*(i0_ac.'*K2*u0_ac+u0_ac.'*K2*i0_ac));
            G_SYNC2=Gsync2/(1-Gsync2*(i0_ac.'*K2*u0_ac+u0_ac.'*K2*i0_ac));
            Gccsc1_iac = 2*(K2*m0_ccsc(1:2)-Gccsc1_f1*K2*i0_cir2)*G_SYNC1*u0_ac.';
            Gccsc2_iac = 2*(K2*m0_ccsc(1:2)-Gccsc2_f1*K2*i0_cir2)*G_SYNC2*u0_ac.';
            %% performing the linear transformation from dq-frame to rotating modified sequence domain
            Gccsc_iac_mpn1 = inv(Tpn2dq) * Tframe_f2 * Gccsc1_iac * inv(Tframe_f1) * Tpn2dq;
            Gccsc_iac_mpn2 = inv(Tpn2dq) * Tframe_f2 * Gccsc2_iac * inv(Tframe_f1) * Tpn2dq;
            %% frequency shifting already included in the frequency-shifted controller
            Gccsc_iac_pn = [Gccsc_iac_mpn1(1,1)     Gccsc_iac_mpn2(1,2)
                            Gccsc_iac_mpn1(2,1)     Gccsc_iac_mpn2(2,2)];
            Gccsc_iac_3 = Tpn2abc * [1 0; 0 0] * Gccsc_iac_pn * Tabc2pn;
            Gccsc_iac3 = Tpn2abc * [0 0; 0 1] * Gccsc_iac_pn * Tabc2pn;
        end
        %% Influence of angle dynamic on CCSC ended
    else
        Gccsc_abc = zeros(3);
        
        Gccsc_iac_abc = zeros(3);
        Gccsc_iac_3 = zeros(3);
        Gccsc_iac3 = zeros(3);
        
        Gwt_uc_sum_z = zeros(3);
        Gwt_uc_diff = zeros(3);
        Gwt_uc_sum_dq = zeros(3);
    end
    
     Hcc1 = Kp_Ireg + Ki_Ireg/(s+j*h*w1-j*w1);
     Hcc2 = Kp_Ireg + Ki_Ireg/(s+j*h*w1+j*w1);
    
    if iMMC.CtrlMode == GFL_APC || iMMC.CtrlMode == GFL_DVC  
        Gcc1 = [ Hcc1       Lff;  -Lff      Hcc1];
        Gcc2 = [ Hcc2       +Lff;  -Lff       Hcc2];
        %H_Pctrl1 = PQCTRL.PI.kp*(1+1/(PQCTRL.PI.Ti*s)); %% here Kp and Ti are given
        H_Pctrl1 = ( Kp_Preg + Ki_Preg/(s+j*h*w1-j*w1) ) * 1; %% for positive sequence
        H_Pctrl2 = ( Kp_Preg + Ki_Preg/(s+j*h*w1+j*w1) ) * 1; %% for negative sequence
        %% Q control with -1 coefficient for the PI control output
        H_Qctrl1 = ( -Kp_Qreg - Ki_Qreg/(s+j*h*w1-j*w1) ) * Q_LPF1 * 1; %% for positive sequence
        H_Qctrl2 = ( -Kp_Qreg - Ki_Qreg/(s+j*h*w1+j*w1) ) * Q_LPF2 * 1; %% for negative sequence

        if iMMC.CtrlMode == GFL_APC
            G_PQi1 = [-H_Pctrl1*u0_ac.'; H_Qctrl1*u0_ac.'*K2];
            G_PQi2 = [-H_Pctrl2*u0_ac.'; H_Qctrl2*u0_ac.'*K2];
            if iMMC.Ctrl.IsVacGFL == 1
                G_PQi1 = [-H_Pctrl1*u0_ac.'; zeros(1,2)];
                G_PQi2 = [-H_Pctrl2*u0_ac.'; zeros(1,2)];
            end
        elseif iMMC.CtrlMode == GFL_DVC
            G_PQi1 = [zeros(1,2); H_Qctrl1*u0_ac.'*K2];
            G_PQi2 = [zeros(1,2); H_Qctrl2*u0_ac.'*K2];
        end
        %Udc0=1;
        Gmac_iac1 = (Hcc1*G_PQi1 - I_LPF1*Gcc1) * Gd * Kpwm/ Udc0;
        Gmac_iac2 = (Hcc2*G_PQi2 - I_LPF2*Gcc2) * Gd * Kpwm/ Udc0;
        %% When neglecting PQ control loops
        % Gmac_iac1 = - I_LPF1*Gcc1 * Gd * Kpwm/ Udc0;
        % Gmac_iac2 = - I_LPF2*Gcc2 * Gd * Kpwm/ Udc0;
        %% performing the linear transformation from dq-frame to rotating modified sequence domain
        Gmac_iac_mpn1 = inv(Tpn2dq) * Gmac_iac1 * Tpn2dq;
        Gmac_iac_mpn2 = inv(Tpn2dq) * Gmac_iac2 * Tpn2dq;
        %% Adding the rotation from SRF to local PoC dq-frame
        if 1
        Tframe_local = [cos(0-theta0) sin(0-theta0); -sin(0-theta0) cos(0-theta0)];
        Gmac_iac_mpn1 = inv(Tpn2dq) * Tframe_local * Gmac_iac1 * inv(Tframe_local) * Tpn2dq;
        Gmac_iac_mpn2 = inv(Tpn2dq) * Tframe_local * Gmac_iac2 * inv(Tframe_local) * Tpn2dq;
        end
    elseif iMMC.CtrlMode == GFM_APC
        %% VSM Modelling
        C2x2=K2;
        %% Active damping denoted by Gad is by default deactivated 
        Gad1=0;% ((s+j*h*w1-j*w1)+VSM.VTFF.OmLead)/((s+j*h*w1-j*w1)+VSM.VTFF.OmLag)*eye(2);
        Gad2=0;% ((s+j*h*w1+j*w1)+VSM.VTFF.OmLead)/((s+j*h*w1+j*w1)+VSM.VTFF.OmLag)*eye(2);
        Grv1=VSM.VC.tau_rv*(s+j*h*w1-j*w1)/(1+VSM.VC.tau_rv*(s+j*h*w1-j*w1)); % mathematical equivalent with approximation
        Grv2=VSM.VC.tau_rv*(s+j*h*w1+j*w1)/(1+VSM.VC.tau_rv*(s+j*h*w1+j*w1));
        Gvc1 = iMMC.Ctrl.IsVIinGFM * (Grv1*VSM.VC.rv*eye(2)-VSM.VC.xv*C2x2);
        Gvc2 = iMMC.Ctrl.IsVIinGFM * (Grv2*VSM.VC.rv*eye(2)-VSM.VC.xv*C2x2);
        %% ADDING RPC
        G_Qe1 = -VSM.ku * 1/(1+VSM.Tu*(s+j*h*w1-j*w1))*Q_LPF1;
        G_Qe2 = -VSM.ku * 1/(1+VSM.Tu*(s+j*h*w1+j*w1))*Q_LPF2;
        Geu1 = [G_Qe1*i0_ac.'*C2x2; zeros(1,2)];
        Geu2 = [G_Qe2*i0_ac.'*C2x2; zeros(1,2)];
        Gei1 = [-G_Qe1*u0_ac.'*C2x2; zeros(1,2)];
        Gei2 = [-G_Qe2*u0_ac.'*C2x2; zeros(1,2)];
        Gad1 = Gad1 + Geu1;
        Gad2 = Gad2 + Geu2;
        Gvc1 = Gvc1 - Gei1;
        Gvc2 = Gvc2 - Gei2;

        LPF_INF = 1;
        %% multiplying LFP and the available power Pm for f-support added in 06.2022
        Gsync1=-1*wb/(s+j*h*w1-j*w1)*(1+VSM.kdf*(s+j*h*w1-j*w1))/(VSM.Ta*(s+j*h*w1-j*w1)+iMMC.Ctrl.IsPfDroopGFM*VSM.kpf*iMMC.Set.Pref)*LPF_INF;
        Gsync2=-1*wb/(s+j*h*w1+j*w1)*(1+VSM.kdf*(s+j*h*w1+j*w1))/(VSM.Ta*(s+j*h*w1+j*w1)+iMMC.Ctrl.IsPfDroopGFM*VSM.kpf*iMMC.Set.Pref)*LPF_INF;
        if 0 %% when the active power is calculated using dq variables
        G_SYNC1=Gsync1/(1-Gsync1*(i0_ac.'*C2x2*u0_ac+u0_ac.'*C2x2*i0_ac));
        G_SYNC2=Gsync2/(1-Gsync2*(i0_ac.'*C2x2*u0_ac+u0_ac.'*C2x2*i0_ac));
        else %% when the active power is calculated using abc variables 
        G_SYNC1=Gsync1;
        G_SYNC2=Gsync2;
        end
        if iMMC.Ctrl.IsCCinGFM == 0 %% without inner loop
            Gmac_iac1 = ( G_SYNC1*(-C2x2*m0_ac/Kpwm*Udc0-Gvc1*I_LPF1*C2x2*i0_ac+Gad1*LPF_INF*C2x2*u0_ac)*u0_ac.'- Gvc1*I_LPF1*eye(2) ) * Gd * Kpwm / Udc0;
            Gmac_iac2 = ( G_SYNC2*(-C2x2*m0_ac/Kpwm*Udc0-Gvc2*I_LPF2*C2x2*i0_ac+Gad2*LPF_INF*C2x2*u0_ac)*u0_ac.'- Gvc2*I_LPF2*eye(2) ) * Gd * Kpwm / Udc0;
        elseif iMMC.Ctrl.IsCCinGFM == 1 %% with inner loop
            %% equivalent to replacing Gvc with Gi2m_tmp and replacing Gad with Gu2m_tmp
            Gu2m_tmp1 = (Kpv_GFM*Hcc1*(Gad1-eye(2)) + eye(2));
            Gu2m_tmp2 = (Kpv_GFM*Hcc2*(Gad2-eye(2)) + eye(2));
            Gi2m_tmp1 = 1*Lff*C2x2 + Kpv_GFM*Hcc1*Gvc1;
            Gi2m_tmp2 = 1*Lff*C2x2 + Kpv_GFM*Hcc2*Gvc2;
            Gmac_iac1 = ( G_SYNC1*(-C2x2*m0_ac/Kpwm*Udc0 -Gi2m_tmp1*I_LPF1*C2x2*i0_ac +Gu2m_tmp1*LPF_INF*C2x2*u0_ac)*u0_ac.' -Gi2m_tmp1*I_LPF1 ) * Gd * Kpwm / Udc0;
            Gmac_iac2 = ( G_SYNC2*(-C2x2*m0_ac/Kpwm*Udc0 -Gi2m_tmp2*I_LPF2*C2x2*i0_ac +Gu2m_tmp2*LPF_INF*C2x2*u0_ac)*u0_ac.' -Gi2m_tmp2*I_LPF2 ) * Gd * Kpwm / Udc0;
        end
        %% DEBUGGING
        % Gmac_iac1 = ( G_SYNC1*(-C2x2*m0_ac)*u0_ac.' ) * Gd ;
        % Gmac_iac2 = ( G_SYNC2*(-C2x2*m0_ac)*u0_ac.' ) * Gd ;
        %% performing the linear transformation from dq-frame to rotating modified sequence domain
        %% and Adding the rotation from SRF to local PoC dq-frame
        d_phi0 = 0 - theta0 * 1;
        Tframe_local = [cos(d_phi0) sin(d_phi0); -sin(d_phi0) cos(d_phi0)];
        Gmac_iac_mpn1 = inv(Tpn2dq) * Tframe_local * Gmac_iac1 * inv(Tframe_local) * Tpn2dq;
        Gmac_iac_mpn2 = inv(Tpn2dq) * Tframe_local * Gmac_iac2 * inv(Tframe_local) * Tpn2dq;
    else
        Gmac_iac_mpn1 = zeros(2);
        Gmac_iac_mpn2 = zeros(2);
    end
    
    %% frequency shifting already included in the frequency-shifted controller
    Gmac_iac_pn = [ Gmac_iac_mpn1(1,1)     0*Gmac_iac_mpn2(1,2)
                    0*Gmac_iac_mpn1(2,1)     Gmac_iac_mpn2(2,2)];
    
    Gmac_iac_P = [      0                   0
                    Gmac_iac_mpn1(2,1)      0   ];
            
    Gmac_iac_N = [  0       Gmac_iac_mpn2(1,2)
                    0               0           ];
    Gmac_iac_abc = Tpn2abc * Gmac_iac_pn * Tabc2pn;
    
    Gmx_out.D = [   Gmac_iac_abc	zeros(3)       zeros(3)         zeros(3)
                    zeros(3)    	Gccsc_abc      Gwt_uc_sum_z*1     Gwt_uc_diff*0    ];
    %% W2 amd W1 are the transfer functions from Uc_sum_dq (n2) and Uc_diff_dq (p1) to energy control output (z0)           
    Gmx_out.W1 = [  zeros(3)	zeros(3)      zeros(3)     zeros(3)
                    zeros(3)    zeros(3)      zeros(3)     Gwt_uc_diff ];
    Gmx_out.W2 = [  zeros(3)	zeros(3)      zeros(3)          zeros(3)
                    zeros(3)    zeros(3)      Gwt_uc_sum_dq     zeros(3)  ];
                
    Gmx_out.L = 1*[     Tpn2abc*Gmac_iac_N*Tabc2pn      zeros(3)       zeros(3)     zeros(3)
                        zeros(3)                        zeros(3)      zeros(3)     zeros(3)    ];
    Gmx_out.R = 1*[     Tpn2abc*Gmac_iac_P*Tabc2pn      zeros(3)       zeros(3)     zeros(3)
                        zeros(3)                        zeros(3)      zeros(3)     zeros(3)    ];
    
    Gmx_out.ccsc3 = [   zeros(3)        zeros(3)       zeros(3)     zeros(3)
                        Gccsc_iac3    	zeros(3)      zeros(3)     zeros(3)    ];
    Gmx_out.ccsc_3 = [  zeros(3)        zeros(3)       zeros(3)     zeros(3)
                        Gccsc_iac_3    	zeros(3)      zeros(3)     zeros(3)    ];
    Gmx_out_vec(k) = Gmx_out;
end

GmxStru.f_3 = Gmx_out_vec(1); GmxStru.f_2 = Gmx_out_vec(2); GmxStru.f_1 = Gmx_out_vec(3); GmxStru.f0 = Gmx_out_vec(4);
GmxStru.f1 = Gmx_out_vec(5); GmxStru.f2 = Gmx_out_vec(6); GmxStru.f3 = Gmx_out_vec(7);