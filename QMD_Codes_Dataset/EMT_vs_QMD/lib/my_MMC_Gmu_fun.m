function GmuStru = my_MMC_Gmu_fun(IN,iMMC)

my_MMC_param

%% Transfer function from u_abc to m_abc

%% m_abc = [mac     mdc]
%% x_abc = [iac_abc   icir*ones(1,3)    uc_abc_cm   uc_abc_dm].' 
%% u_abc = [uac_abc     udc*ones(1,3)].'

w1 = 2 * pi * 50;

H_offset=1e-4; %% pu of harmonic order

u0_ac = IN.u0_ac;
i0_ac = IN.i0_ac;
m0_ac = IN.m0_ac;
m0_ccsc = IN.m0_cir;
i0_cir2 = IN.i0_cir2;
theta0 = IN.theta0;
U1d=u0_ac(1); U1q=u0_ac(2); %%% Us_vec=U1d+jU1q, measured values
I1d=i0_ac(1); I1q=i0_ac(2); %%% Ig1_vec=Ig1d+jIg1q, measured values
M1d=m0_ac(1); M1q=m0_ac(2);

M2d=m0_ccsc(1); M2q=m0_ccsc(2);
I2d_cir=i0_cir2(1); I2q_cir=i0_cir2(2); %% approximate 0
% U0=sqrt(U1d^2+U1q^2);
% P0 = i0_ac.'*u0_ac;

K2 = [0   1;    -1   0];

h_vec = [-3 -2 -1 0 1 2 3];

for k=1:1:length(h_vec)
    h = h_vec(k);
    s = IN.s;
    if s-1j*w1+j*h*w1 == 0 || s+1j*w1+j*h*w1 == 0
        s = s-1j*H_offset*w1;
    end
    
    %format long
    % Gd = 1/(1+1.5*s*Ts);%% inaccurate by large delay
    % Gd = exp(-1.5*s*Ts);
    % Gd1 = 1/(1+1.5*(s+j*h*w1-j*w1)*Ts);
    % Gd2 = 1/(1+1.5*(s+j*h*w1+j*w1)*Ts);
    Gd = exp(-1.0*s*Ts)*(1-exp(-1.5*s*Ts))/(1.5*s*Ts);
    
    H_LPF1 = (2*pi*Fn_filter)^2 / ( (s+j*h*w1-j*w1)^2 + 2*Zeta_filter*(2*pi*Fn_filter)*(s+j*h*w1-j*w1) + (2*pi*Fn_filter)^2 );
    H_LPF2 = (2*pi*Fn_filter)^2 / ( (s+j*h*w1+j*w1)^2 + 2*Zeta_filter*(2*pi*Fn_filter)*(s+j*h*w1+j*w1) + (2*pi*Fn_filter)^2 );
    % C_LPF1 = (2*pi*Fn_filter)^2 / ( (s+j*h*w1-(-2j*w1))^2 + 2*Zeta_filter*(2*pi*Fn_filter)*(s+j*h*w1-(-2j*w1)) + (2*pi*Fn_filter)^2 );
    % C_LPF2 = (2*pi*Fn_filter)^2 / ( (s+j*h*w1+(-2j*w1))^2 + 2*Zeta_filter*(2*pi*Fn_filter)*(s+j*h*w1+(-2j*w1)) + (2*pi*Fn_filter)^2 );

    I_LPF1 = 1; I_LPF2 = 1;
    C_LPF1 = 1; C_LPF2 = 1;

    Q_LPF1 = (2*pi*Fn_filter/10)^2 / ( (s+j*h*w1-j*w1)^2 + 2*Zeta_filter*(2*pi*Fn_filter/10)*(s+j*h*w1-j*w1) + (2*pi*Fn_filter/10)^2 );
    Q_LPF2 = (2*pi*Fn_filter/10)^2 / ( (s+j*h*w1+j*w1)^2 + 2*Zeta_filter*(2*pi*Fn_filter/10)*(s+j*h*w1+j*w1) + (2*pi*Fn_filter/10)^2 );
    %% Used only in the CC loop as FF
    U_LPF1 = 1 / (1 + 1/(2*pi*Fn_filter/100)*(s+j*h*w1-j*w1));
    U_LPF2 = 1 / (1 + 1/(2*pi*Fn_filter/100)*(s+j*h*w1+j*w1));
    
    %% frequency-shifted CCSC transfer function used for modelling VSM angle perturbations
    Hccsc1 = Kp_Iccsc + Ki_Iccsc/(s+j*h*w1-(-2j*w1)); %% for positive sequence
    Hccsc2 = Kp_Iccsc + Ki_Iccsc/(s+j*h*w1+(-2j*w1)); %% por negative sequence
    %% non-diagonal terms of the following matrix is obtained using the p/n relationships 
    Gccsc1 = [  Hccsc1      -2*Lff
                +2*Lff      Hccsc1  ] * C_LPF1 * Gd / Udc0;
    %% The following transfer function is formulated based on the CCSC control structure        
    Gccsc2 = [  Hccsc2      -2*Lff
                2*Lff      Hccsc2  ] * C_LPF2 * Gd / Udc0;
    
    Gmac_vdc1 = zeros(2); Gmac_vdc2 = zeros(2); gmac_vdc = zeros(2,1);
    
    Hcc1 = Kp_Ireg + Ki_Ireg/(s+j*h*w1-j*w1);
    Hcc2 = Kp_Ireg + Ki_Ireg/(s+j*h*w1+j*w1);

    %% Main Controller
    if iMMC.CtrlMode == GFM_VF
        %% frequency-shifted AC voltage control from uac_abc to mac
        Hvac1 = Kp_VACreg + Ki_VACreg/(s+j*h*w1-j*w1); %% for positive sequence
        Hvac2 = Kp_VACreg + Ki_VACreg/(s+j*h*w1+j*w1); %% for negative sequence
        Kdd = Lff * 0;
        %% transfer frunction between mac_dq output and and Uac_dq input
        Gmac_vac1 = -[Hvac1     +Kdd
                      -Kdd       Hvac1]  * Gd / Udc0; %* H_LPF1
        Gmac_vac2 = -[Hvac2     +Kdd
                      -Kdd       Hvac2]  * Gd / Udc0; %* H_LPF2
        Gccsc_abc = zeros(3);
        Gccsc_abc_3 = zeros(3);
        Gccsc_abc3 = zeros(3);
    elseif iMMC.CtrlMode == GFL_APC || iMMC.CtrlMode == GFL_DVC
        %% Applies to PLL-synchronized control modes
        Hpll1 = Kp_PLL + Ki_PLL/(s+j*h*w1-j*w1); % positive sequence;
        Hpll2 = Kp_PLL + Ki_PLL/(s+j*h*w1+j*w1); % positive sequence;
        Gpll1 = Hpll1/((s+j*h*w1-j*w1)+Hpll1*U1d);
        Gpll2 = Hpll2/((s+j*h*w1+j*w1)+Hpll2*U1d);
        Gpllu1=[1  Gpll1*U1q; 0  1-Gpll1*U1d];
        Gpllu2=[1  Gpll2*U1q; 0  1-Gpll2*U1d];
        Gplli1=[0  Gpll1*I1q; 0  -Gpll1*I1d];
        Gplli2=[0  Gpll2*I1q; 0  -Gpll2*I1d];
        Gpllm1=[0  -Gpll1*M1q; 0  Gpll1*M1d];
        Gpllm2=[0  -Gpll2*M1q; 0  Gpll2*M1d];
        
        Hcc = Kp_Ireg + Ki_Ireg/(s+j*h*w1);
        
        Gcc1 = [ Hcc1       Lff;  -Lff      Hcc1];
        Gcc2 = [ Hcc2       +Lff;  -Lff       Hcc2];
        %H_Pctrl1 = PQCTRL.PI.kp*(1+1/(PQCTRL.PI.Ti*s)); %% here Kp and Ti are given
        H_Pctrl1 = ( Kp_Preg + Ki_Preg/(s+j*h*w1-j*w1) ) * 1; %% for positive sequence
        H_Pctrl2 = ( Kp_Preg + Ki_Preg/(s+j*h*w1+j*w1) ) * 1; %% for negative sequence
        %% Q control with -1 coefficient for the PI control output
        H_Qctrl1 = ( -Kp_Qreg - Ki_Qreg/(s+j*h*w1-j*w1) ) * Q_LPF1 * 1; %% for positive sequence
        H_Qctrl2 = ( -Kp_Qreg - Ki_Qreg/(s+j*h*w1+j*w1) ) * Q_LPF2 * 1; %% for negative sequence
        %% For PV mode: LPF added to the Vac Controller
        if iMMC.Ctrl.IsVacGFL == 1
            GFL_H_Vac1 = ( Kp_VAC_GFL + Ki_VAC_GFL/(s+j*h*w1-j*w1) ) * U_LPF1; %% for positive sequence
            GFL_H_Vac2 = ( Kp_VAC_GFL + Ki_VAC_GFL/(s+j*h*w1+j*w1) ) * U_LPF2; %% for negative sequence
        end

        if iMMC.CtrlMode == GFL_APC
            G_PQu1 = [-H_Pctrl1*i0_ac.'; -H_Qctrl1*i0_ac.'*K2];
            G_PQu2 = [-H_Pctrl2*i0_ac.'; -H_Qctrl2*i0_ac.'*K2];
            G_PQi1 = [-H_Pctrl1*u0_ac.'; H_Qctrl1*u0_ac.'*K2];
            G_PQi2 = [-H_Pctrl2*u0_ac.'; H_Qctrl2*u0_ac.'*K2];
            if iMMC.Ctrl.IsVacGFL == 1
                G_PQu1 = [-H_Pctrl1*i0_ac.'; GFL_H_Vac1*u0_ac.'/sqrt(u0_ac(1)^2+u0_ac(2)^2)];
                G_PQu2 = [-H_Pctrl2*i0_ac.'; GFL_H_Vac2*u0_ac.'/sqrt(u0_ac(1)^2+u0_ac(2)^2)];
                G_PQi1 = [-H_Pctrl1*u0_ac.'; zeros(1,2)];
                G_PQi2 = [-H_Pctrl2*u0_ac.'; zeros(1,2)];
            end
        elseif iMMC.CtrlMode == GFL_DVC
            G_PQu1 = [zeros(1,2); -H_Qctrl1*i0_ac.'*K2];
            G_PQu2 = [zeros(1,2); -H_Qctrl2*i0_ac.'*K2];
            G_PQi1 = [zeros(1,2); H_Qctrl1*u0_ac.'*K2];
            G_PQi2 = [zeros(1,2); H_Qctrl2*u0_ac.'*K2];
            %% Vdc control transfer function: refer to "Impedance Modeling and Analysis of MMC-HVDC for OWF Integration"
            Hvdc = Kp_VDCreg + Ki_VDCreg/(s+j*h*w1);
            %% p.u. udc should be carefully dealt with
            Gvdc_dq = Ubase/Vnom_dc * Hcc * Hvdc * Gd * Kpwm / Udc0;
            gmac_vdc = Gvdc_dq * [1 0].';
        end
        Gmac_vac1 = ( (Hcc1*G_PQu1+U_LPF1*eye(2))*Gpllu1 + Gpllm1/Kpwm*Udc0 + (Hcc1*G_PQi1-I_LPF1*Gcc1)*Gplli1 ) * Gd * Kpwm / Udc0;
        Gmac_vac2 = ( (Hcc2*G_PQu2+U_LPF2*eye(2))*Gpllu2 + Gpllm2/Kpwm*Udc0 + (Hcc2*G_PQi2-I_LPF2*Gcc2)*Gplli2 ) * Gd * Kpwm / Udc0;
        %% When neglecting PQ control loops
        % Gmac_vac1 = ( 1*U_LPF1*Gpllu1 + 1*Gpllm1*Udc0/Kpwm - 1*I_LPF1*Gcc1*Gplli1 ) * Gd * Kpwm / Udc0;
        % Gmac_vac2 = ( 1*U_LPF2*Gpllu2 + 1*Gpllm2*Udc0/Kpwm - 1*I_LPF2*Gcc2*Gplli2 ) * Gd * Kpwm / Udc0;
        %% Also neglecting PLL OTHER THAN outer-loop dynamics 
        % Gmac_vac1 = eye(2) * U_LPF1 * Gd1 * Kpwm / Udc0;
        % Gmac_vac2 = eye(2) * U_LPF2 * Gd2 * Kpwm / Udc0;

        %% CCSC response to PLL angle perturbation is negligible, as the steady-state icir is very small after implementing CCSC control 
        if iMMC.Ctrl.CCSC.Flag == 1
            %% Add PLL angle dynamics
            Gdq2_plli1=[0  -2*Gpll1*I2q_cir; 0  2*Gpll1*I2d_cir];
            Gdq2_plli2=[0  -2*Gpll2*I2q_cir; 0  2*Gpll2*I2d_cir];
            Gdq2_pllm1=[0  -2*Gpll1*M2q; 0  2*Gpll1*M2d];
            Gdq2_pllm2=[0  -2*Gpll2*M2q; 0  2*Gpll2*M2d];
            Gccsc1 = Gccsc1 * Gdq2_plli1 - Gdq2_pllm1;
            Gccsc2 = Gccsc2 * Gdq2_plli2 - Gdq2_pllm2;
            %% performing the linear transformation from dq-frame to rotating modified sequence domain
            Gccsc_mpn1 = inv(Tpn2dq) * Gccsc1 * Tpn2dq;
            Gccsc_mpn2 = inv(Tpn2dq) * Gccsc2 * Tpn2dq;
                if 1 %% adding the following rotation will not change the results
                theta0_f2 = -theta0*2;
                Tframe_local = [cos(0-theta0_f2) sin(0-theta0_f2); -sin(0-theta0_f2) cos(0-theta0_f2)];
                Gccsc_mpn1 = inv(Tpn2dq) * Tframe_local * Gccsc1 * inv(Tframe_local) * Tpn2dq;
                Gccsc_mpn2 = inv(Tpn2dq) * Tframe_local * Gccsc2 * inv(Tframe_local) * Tpn2dq;
                end
            %% frequency shifting already included in the frequency-shifted controller
            Gccsc_pn = [Gccsc_mpn1(1,1)     Gccsc_mpn2(1,2)
                        Gccsc_mpn1(2,1)     Gccsc_mpn2(2,2)];
            Gccsc_abc = Tpn2abc * Gccsc_pn * Tabc2pn;
            %% following 2 lines modified from 0 on 30.01.2024
            Gccsc_abc_3 = Tpn2abc * [1 0; 0 0] * Gccsc_pn * Tabc2pn; %% positive sequence
            Gccsc_abc3 = Tpn2abc * [0 0; 0 1] * Gccsc_pn * Tabc2pn; %% negative sequence
        else
            Gccsc_abc = zeros(3);
            Gccsc_abc_3 = zeros(3);
            Gccsc_abc3 = zeros(3);
        end
    elseif iMMC.CtrlMode == GFM_APC
        %% VSM Modelling
        C2x2 = K2;
        %% Active damping denoted by Gad is by default deactivated 
        Gad1=0; %((s+j*h*w1-j*w1)+VSM.VTFF.OmLead)/((s+j*h*w1-j*w1)+VSM.VTFF.OmLag)*eye(2);
        Gad2=0; %((s+j*h*w1+j*w1)+VSM.VTFF.OmLead)/((s+j*h*w1+j*w1)+VSM.VTFF.OmLag)*eye(2);
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
        %% multiplying LFP and multiplying available Pm for f-support in 06.2022
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
            Gmac_vac1 = ( G_SYNC1*(-C2x2*m0_ac/Kpwm*Udc0 -Gvc1*I_LPF1*C2x2*i0_ac+Gad1*LPF_INF*C2x2*u0_ac)*i0_ac.'+ Gad1*LPF_INF*eye(2) ) * Gd * Kpwm / Udc0;
            Gmac_vac2 = ( G_SYNC2*(-C2x2*m0_ac/Kpwm*Udc0 -Gvc2*I_LPF2*C2x2*i0_ac+Gad2*LPF_INF*C2x2*u0_ac)*i0_ac.'+ Gad2*LPF_INF*eye(2) ) * Gd * Kpwm / Udc0;
        elseif iMMC.Ctrl.IsCCinGFM == 1 %% with inner loop
            %% equivalent to replacing Gvc with Gi2m_tmp and replacing Gad with Gu2m_tmp
            Gu2m_tmp1=(Kpv_GFM*Hcc1*(Gad1-eye(2))+eye(2));
            Gu2m_tmp2=(Kpv_GFM*Hcc2*(Gad2-eye(2))+eye(2));
            Gi2m_tmp1=Lff*C2x2+Kpv_GFM*Hcc1*Gvc1;
            Gi2m_tmp2=Lff*C2x2+Kpv_GFM*Hcc2*Gvc2;
            Gmac_vac1 = ( G_SYNC1*(-C2x2*m0_ac/Kpwm*Udc0 -Gi2m_tmp1*I_LPF1*C2x2*i0_ac +Gu2m_tmp1*LPF_INF*C2x2*u0_ac)*i0_ac.'+ Gu2m_tmp1*LPF_INF ) * Gd * Kpwm / Udc0;
            Gmac_vac2 = ( G_SYNC2*(-C2x2*m0_ac/Kpwm*Udc0 -Gi2m_tmp2*I_LPF2*C2x2*i0_ac +Gu2m_tmp2*LPF_INF*C2x2*u0_ac)*i0_ac.'+ Gu2m_tmp2*LPF_INF ) * Gd * Kpwm / Udc0;
        end
        % DEBUGGING
        % Gmac_vac1 = ( G_SYNC1*(-C2x2*m0_ac )*i0_ac.' ) * Gd ;
        % Gmac_vac2 = ( G_SYNC2*(-C2x2*m0_ac )*i0_ac.' ) * Gd ;
    
        %% CCSC response to VSM angle perturbation is negligible ???
        if iMMC.Ctrl.CCSC.Flag == 1
            %% Adding VSM angle dynamics
            %% frequency-shifted CCSC transfer function used for modelling VSM angle perturbations
            Hccsc1_f1 = Kp_Ireg + Ki_Ireg/(s+j*h*w1-(1j*w1)); %% for positive sequence
            Hccsc2_f1 = Kp_Ireg + Ki_Ireg/(s+j*h*w1+(1j*w1)); %% por negative sequence
            % %% FOR DEBUGGING
            % Hccsc = Kp_Ireg + Ki_Ireg/(s+j*h*w1); %% for positive sequence
            % Hccsc1_f1 = Hccsc; Hccsc2_f1 = Hccsc;
            %% non-diagonal terms of the following matrix is obtained using the p/n relationships 
            Gccsc1_f1 = [  Hccsc1_f1      -2*Lff
                        +2*Lff      Hccsc1_f1  ] * C_LPF1 * Gd / Udc0;
            %% The following transfer function is formulated based on the CCSC control structure        
            Gccsc2_f1 = [  Hccsc2_f1      -2*Lff
                        2*Lff      Hccsc2_f1  ] * C_LPF2 * Gd / Udc0;
            Gccsc1_f1 = 2*(K2*m0_ccsc(1:2)-Gccsc1_f1*K2*i0_cir2)*G_SYNC1*i0_ac.';
            Gccsc2_f1 = 2*(K2*m0_ccsc(1:2)-Gccsc2_f1*K2*i0_cir2)*G_SYNC2*i0_ac.';
            %% performing the linear transformation from dq-frame to rotating modified sequence domain
            Gccsc_mpn1 = inv(Tpn2dq) * Gccsc1_f1 * Tpn2dq;
            Gccsc_mpn2 = inv(Tpn2dq) * Gccsc2_f1 * Tpn2dq;
            if 1 %% adding the following rotation MAY not change the results
                theta0_f2 = -theta0*2;
                Tframe_local = [cos(0-theta0_f2) sin(0-theta0_f2); -sin(0-theta0_f2) cos(0-theta0_f2)];
                Gccsc_mpn1 = inv(Tpn2dq) * Tframe_local * Gccsc1_f1 * inv(Tframe_local) * Tpn2dq;
                Gccsc_mpn2 = inv(Tpn2dq) * Tframe_local * Gccsc2_f1 * inv(Tframe_local) * Tpn2dq;
            end
            %% frequency shifting already included in the frequency-shifted controller
            Gccsc_pn = [Gccsc_mpn1(1,1)     Gccsc_mpn2(1,2)
                        Gccsc_mpn1(2,1)     Gccsc_mpn2(2,2)];
            Gccsc_abc = Tpn2abc * Gccsc_pn * Tabc2pn;
            Gccsc_abc_3 = Tpn2abc * [1 0; 0 0] * Gccsc_pn * Tabc2pn; %% positive sequence
            Gccsc_abc3 = Tpn2abc * [0 0; 0 1] * Gccsc_pn * Tabc2pn; %% negative sequence
        else
            Gccsc_abc = zeros(3);
            Gccsc_abc_3 = zeros(3);
            Gccsc_abc3 = zeros(3);
        end
    end
    
    %% Adding the rotation from SRF to local system/PoC dq-frame
    d_phi0 = 0 - theta0 * 1;
    Tframe_local = [cos(d_phi0) sin(d_phi0); -sin(d_phi0) cos(d_phi0)];
    Gmac_vac_mpn1 = inv(Tpn2dq) * (Tframe_local * Gmac_vac1 * inv(Tframe_local)) * Tpn2dq;
    Gmac_vac_mpn2 = inv(Tpn2dq) * (Tframe_local * Gmac_vac2 * inv(Tframe_local)) * Tpn2dq;
    
    %% udc control output rotation -- Only the ac varibles should be rotated to the global dq-frame
    gmac_vdc = Tdq2pn * Tframe_local * gmac_vdc;
    
    %% frequency shifting already included in the frequency-shifted controller
    Gmac_vac_pn = [ Gmac_vac_mpn1(1,1)     0*Gmac_vac_mpn2(1,2)
                    0*Gmac_vac_mpn1(2,1)     Gmac_vac_mpn2(2,2)];
            
    Gmac_vac_P = [      0                   0
                    Gmac_vac_mpn1(2,1)      0   ];
            
    Gmac_vac_N = [  0       Gmac_vac_mpn2(1,2)
                    0               0           ];
            
    Gmac_vac_abc = Tpn2abc * Gmac_vac_pn * Tabc2pn;
    
    if iMMC.CtrlMode == GFL_DVC
    %% udc: Note that Tdq2pn*[1; 0] = [0.5 0.5], if neglecting the steady-state angle rotation, Gmac_p/n = 0.5 * Gvdc
        Gmac_vdc_abc1 = diag( Tpn2abc * ([1 0; 0 0] * gmac_vdc) );
        Gmac_vdc_abc_1 = diag( Tpn2abc * ([0 0; 0 1] * gmac_vdc) );
        % Gmac_vdc_abc1 = diag( Tp_abc * gmac_vdc(1) );
        % Gmac_vdc_abc_1 = diag( Tn_abc * gmac_vdc(2) );
    else
        Gmac_vdc_abc1 = zeros(3);
        Gmac_vdc_abc_1 = zeros(3);
    end
    
    %% if u = [vac vdc*ones(1,3)]
    %% Output diagonal,left and right terms
    % -----------------------------------------------------------
    Gmu_out.D = [   Gmac_vac_abc        zeros(3)
                    zeros(3)            zeros(3)    ];
     
    Gmu_out.dc1 = [  zeros(3)        Gmac_vdc_abc1
                    zeros(3)            zeros(3)    ];
    Gmu_out.dc_1 = [  zeros(3)        Gmac_vdc_abc_1
                    zeros(3)            zeros(3)    ];
    %
    Gmu_out.L = 1*[	Tpn2abc*Gmac_vac_N*Tabc2pn      zeros(3)
                    zeros(3)                        zeros(3) ];
    Gmu_out.R = 1*[ Tpn2abc*Gmac_vac_P*Tabc2pn      zeros(3)
                    zeros(3)                        zeros(3) ];
    % -----------------------------------------------------------
    Gmu_out.ccsc3 = [   zeros(3)        zeros(3)
                        Gccsc_abc3      zeros(3)    ];
    Gmu_out.ccsc_3 = [	zeros(3)        zeros(3)
                        Gccsc_abc_3     zeros(3)    ];
    Gmu_out_vec(k) = Gmu_out;
end

GmuStru.f_3 = Gmu_out_vec(1); GmuStru.f_2 = Gmu_out_vec(2); GmuStru.f_1 = Gmu_out_vec(3); GmuStru.f0 = Gmu_out_vec(4);
GmuStru.f1 = Gmu_out_vec(5); GmuStru.f2 = Gmu_out_vec(6); GmuStru.f3 = Gmu_out_vec(7);