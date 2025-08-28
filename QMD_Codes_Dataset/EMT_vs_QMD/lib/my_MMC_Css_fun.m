function C_out = my_MMC_Css_fun(x_SS)

my_MMC_param


for k=1:1:7
    h = k-4;
    %% x_abc = [iac_abc   icir*ones(1,3)    uc_abc_cm   uc_abc_dm].' 
    index = 3*12+h*12;
    Iac0_abc = x_SS( (index+1):(index+3) );
    Icir0_abc = x_SS( (index+4):(index+6) );
    Uc0_sum = x_SS( (index+7):(index+9) );
    Uc0_diff = x_SS( (index+10):(index+12) );
    
    
    if IS_VNO_INCLUDED == 1
        Uc0_sum_M = [2/3 -1/3 -1/3; -1/3 2/3 -1/3; -1/3 -1/3 2/3]*diag(Uc0_sum);
        Uc0_diff_M = -[-1/3 1/6 1/6; 1/6 -1/3 1/6; 1/6 1/6 -1/3]*diag(Uc0_diff)*2;
    
        Ctemp =  [  Uc0_sum_M/(Larm_pu/wb)                  -Uc0_diff_M/2/(Larm_pu/wb)
                    Uc0_diff/2/(Larm_pu/wb).*eye(3)         -Uc0_sum/4/(Larm_pu/wb).*eye(3)
                    -Nb_PM/(Csm_pu/wb)*Iac0_abc.*eye(3)     Nb_PM/(Csm_pu/wb)*Icir0_abc.*eye(3)
                    -2*Nb_PM/(Csm_pu/wb)*Icir0_abc.*eye(3)  Nb_PM/2/(Csm_pu/wb)*Iac0_abc.*eye(3) ];
    else
        Ctemp =  [  Uc0_sum/(Larm_pu/wb).*eye(3)            -Uc0_diff/2/(Larm_pu/wb).*eye(3)
                    Uc0_diff/2/(Larm_pu/wb).*eye(3)         -Uc0_sum/4/(Larm_pu/wb).*eye(3)
                    -Nb_PM/(Csm_pu/wb)*Iac0_abc.*eye(3)     Nb_PM/(Csm_pu/wb)*Icir0_abc.*eye(3)
                    -2*Nb_PM/(Csm_pu/wb)*Icir0_abc.*eye(3)  Nb_PM/2/(Csm_pu/wb)*Iac0_abc.*eye(3) ];
    end

    switch h
        case -3
            C_out.C_3 = Ctemp;
        case -2
            C_out.C_2 = Ctemp;
        case -1
            C_out.C_1 = Ctemp;
        case 0
            C_out.C0 = Ctemp;
        case 1
            C_out.C1 = Ctemp;
        case 2
            C_out.C2 = Ctemp;
        case 3
            C_out.C3 = Ctemp;
        otherwise
            disp('h value beyond [-3 3]')
    end
    %% switch ended here
end
