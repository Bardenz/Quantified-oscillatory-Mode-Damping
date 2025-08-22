%% Notes
% modal-switch correction technique to improve: e.g. the 1922Hz mode is lost is in case118

function RMA_out = my_RMA_fun(f_buf,Ykk_buf, RMA_PLOT_SIGN)
%format shortG
warning('off','all')

if exist("RMA_PLOT_SIGN","var") == 0
    RMA_PLOT_SIGN = 1;
end

IF_PLOT_PF_FLAG = 0;
IF_PLOT_ZM_FLAG = RMA_PLOT_SIGN;
IF_PLOT_PEAK_FLAG = 0;
IF_INCLUDE_ENVELOP = 1;

ZM_LOW_LIMIT = 1;% in abs value? How to generalize it?
DAMPING_UP_LIMIT = 0.25; %% high limit may wrongly treat minor variations as modes

%% --------------- RMA Using Frequency-Sweep Model ---------------------------------
[~, vec_y] = size(Ykk_buf);
nbus = sqrt(vec_y);
display(['nbus = ' num2str(nbus)]);
%% set mStart_NO=nbus+1 to get Envelop, otherwise set mStart_NO=1
mStart_NO = 0*nbus+1;
for k=1:1:length(f_buf)

    s = j*2*pi*f_buf(k);

    % if f_buf(k)-50 == 0
    %     s = j*2*pi*f_buf(k-1);
    % end

    Ykk = reshape(Ykk_buf(k,:),[nbus,nbus]).';

    %[nbus,xxx] = size(Ykk);
    % Ykk_array = reshape(Ykk.',1,[]);
    % Ykk_buf(k,:)= Ykk_array;
    
    %% MAKE SURE the nodal and sequence singals of u,i and um,im for i*Ykk=u, im*Ym=um are in the same locations 
    [T,Ym,L]=eig(Ykk);
    ym=diag(Ym);
    zm=(1./diag(Ym)).';
    
    Ym_buf(k,:)=ym;
    Zm_buf(k,:)=zm;
    %% generating the envelop of Zm frequency-spectrum
    [~,Zm_max_index] = max(abs(Zm_buf(k,:)));
    Zm_max = Zm_buf(k,Zm_max_index);
    
    mNum=nbus;
    Zkk=inv(Ykk);
    Zk_buf(k,:)=diag(Zkk)';%%%buffering driving point impedance
    
    %% CORRECTING MODAL SWITCH CAUSED BY USING eig() TO CALCULATE Ym
    %% improper setting may cause the missing of curve peaks, e.g. when CRT_FACTOR > 0.8
    CRT_FACTOR = 0.95; % default 0.95
    if k==1
        Tbuf=T;
    elseif k>1
        P=T'*Tbuf; Tbuf=T; Tbuf1=T;
        for column=1:1:mNum
            if abs(P(column,column))>CRT_FACTOR
                continue;
            end
            for row=1:1:mNum
                if (abs(P(row,column))>CRT_FACTOR) && (row ~= column)
                    Zm_buf(k,column)=zm(row);
                    Ym_buf(k,column)=ym(row);
                    Tbuf1(:,column)=Tbuf(:,row);
                    break;
                end
            end
        end
        Tbuf=Tbuf1;
    end
    %% ENDING FOR EACH SEQUENCE COMPONENT
    Zm_max_buf(k) = Zm_max;
end
%% Adding the envelop line for checking if some curve peaks are missed
Zm_buf(:,nbus+1)=Zm_max_buf;
Ym_buf(:,nbus+1)=1./Zm_max_buf;

if IF_PLOT_ZM_FLAG == 1

%figure,
%% DISPLAYING MODAL IMPEDANCE CURVE
subplot(3,1,1)
%line_style=['r' 'b' 'm' 'g' 'c' 'k' 'r' 'b' 'm' 'g' 'c' 'k' 'r' 'b' 'm' 'g' 'c' 'k' ];
for m=1:mNum
    line_style(m)="-";
end
line_style(m+1:m+3)="r-";

for m=mStart_NO:mNum+IF_INCLUDE_ENVELOP
    %plot(f_buf,mag2db(abs(Zm_buf(:,m))),line_style(m)),hold on,grid on,
    semilogx(f_buf,(abs(Zm_buf(:,m))),line_style(m)),hold on,grid on,
    %plot(f_buf,(abs(Zm_buf(:,m))),line_style(m)),hold on,grid on,
end
%xlabel('{\it f} in Hz'),
%ylabel('|{\it Z}_i| in pu'),
ylabel('Mag in Ohm'),
%ylabel('Mag in Ohm'),
%set(gca,'xticklabel',[]),
title('Modal Impedances'),% \itZ_i (i=1,2,3)'),
%legend('|{\itZ}_m|'); legend boxoff;

subplot(3,1,2)
for m=mStart_NO:mNum+IF_INCLUDE_ENVELOP
    semilogx(f_buf,real(Zm_buf(:,m)),line_style(m)), hold on,grid on,
    %semilogx(f_buf,angle(Zm_buf(:,m))*180/pi,line_style(m)), hold on,grid on,
    %plot(f_buf,real(Zm_buf(:,m)),line_style(m)), hold on,grid on,
end
ylabel(['Real in Ohm']);
%ylabel(['Phase in deg.']);
%ylabel('{\it R}_m & {\it X}_m in Ohm'),
%xlabel('{\it f} in Hz'),
%set(gca,'xticklabel',[]),
%legend('{\it R}_m_1','{\it X}_m_1'); legend boxoff;

subplot(3,1,3)
for m=mStart_NO:mNum+IF_INCLUDE_ENVELOP
    semilogx(f_buf,imag(Zm_buf(:,m)),line_style(m)), hold on,grid on,
    %plot(f_buf,imag(Zm_buf(:,m)),line_style(m)), hold on,grid on,
end
%ylabel(['\itRe (\itZ\rm_i) in pu']);
ylabel(['Imag. in Ohm']);
%ylabel('{\it R}_m & {\it X}_m in Ohm'),
xlabel('{\it f} in Hz'),
%legend('{\it R}_m_1','{\it X}_m_1'); legend boxoff;
%legend boxoff;
if 0
xticks([0.1 1 10 100])
xticklabels({'0.1','1','10','100'}),xlim([0.1 200]),
end

end

%% Identifying resonance points
res_no=0; Mode_buf=[];
modal_num=nbus;
%for k=2:1:Hmax/Hstep-10
length_Ykk = length(Ykk_buf(:,1));
for k=2:1:length_Ykk-1
    for m=mStart_NO:1:modal_num+1
        if  ( abs(Zm_buf(k,m)) - abs(Zm_buf(k-1,m)) > 1e-12 ) &&  ( abs(Zm_buf(k,m))- abs(Zm_buf(k+1,m)) > 1e-12 )
        %if   ( imag(Zm_buf(k-1,m)) * imag(Zm_buf(k,m)) <= 0 && imag(Zm_buf(k-1,m))-imag(Zm_buf(k,m))~=0) % imaginary part zero-crossing
        %if  ( abs(imag(Ym_buf(k-1,m))) - abs(imag(Ym_buf(k,m))) > 0e-9 ) &&  ( abs(imag(Ym_buf(k+1,m)))- abs(imag(Ym_buf(k,m))) > 0e-9 )
        %if  ( abs(Ym_buf(k-1,m)) - abs(Ym_buf(k,m)) > 1e-9 ) &&  ( abs(Ym_buf(k+1,m))- abs(Ym_buf(k,m)) > 1e-9 ) 
            if abs(Zm_buf(k,m)) < ZM_LOW_LIMIT
                continue;
            end
            %% removing the modes that not appear at the envelop curve and satisfying larger than 5*ZM_LOW_LIMIT
            if abs(Zm_buf(k,m)) < abs(Zm_buf(k,nbus+1)) && real(Zm_buf(k,m)) < real(Zm_buf(k,nbus+1)) ...
                    && imag(Zm_buf(k,m)) < imag(Zm_buf(k,nbus+1)) && real(Zm_buf(k,m))>0 && abs(Zm_buf(k,m)) < 5 * ZM_LOW_LIMIT
                continue;
            end
            [mode_num,~] = size(Mode_buf);
            if m==modal_num+1 && mode_num ~= 0
                if_mode_exist = find(Mode_buf(:,2)==Zm_buf(k,m),1);
                %disp(['If the mode already exist: ' num2str(if_mode_exist)]);
                if if_mode_exist ~= 0
                    continue;
                else
                    disp(['-------------------------------------------------------']);
                    disp(['Missing mode found in envelop curve: ' num2str(f_buf(k)) ' Hz']);
                    disp(['-------------------------------------------------------']);
                end
            end

            %% adaptive bandwidth kBW corresponding to full-width at kBW=0.707 maximum, i.e. -3dB magnitude decrease in dB
            kBW_max = 0.95; kBW_step = 0.05;
            kBW = kBW_max; %% default is 1/sqrt(2)
            MAX_ZETA = DAMPING_UP_LIMIT;
            low_index=0; up_index=0; low_index_kBWmax=0; up_index_kBWmax=0;
            for cur_no=1:1:round(kBW_max/kBW_step)
                % get the min. frequency to consider: where mag-curve rises
                for k2=max(round((1-MAX_ZETA)*k),2):1:k
                    if abs(Zm_buf(k2,m)) - abs(Zm_buf(k,m))*kBW < 1e-12 &&  abs(Zm_buf(k2+1,m))- abs(Zm_buf(k,m))*kBW > 1e-12
                        low_index = k2;
                        %disp(['kBW: ' num2str(kBW) ', f_res-: ' num2str(f_buf(low_index)) ', f_res: ' num2str(f_buf(k)) ' Hz']);
                        if kBW == kBW_max
                            low_index_kBWmax = k2;
                        end
                        break;
                    end
                end
                %% considering very small curve-drop in one side 
                for k2=k:1:min(round((1+MAX_ZETA)*k),length_Ykk-1)
                    if abs(Zm_buf(k2,m)) - abs(Zm_buf(k,m))*kBW > 1e-12 && abs(Zm_buf(k2+1,m))- abs(Zm_buf(k,m))*kBW < 1e-12
                        up_index = k2;
                        %disp(['kBW: ' num2str(kBW) ', f_res+: ' num2str(f_buf(up_index)) ' Hz, f_res: ' num2str(f_buf(k))]);
                        if kBW == kBW_max
                            up_index_kBWmax = k2;
                        end
                        break;
                    end
                end
    
                if (low_index ~= 0 && low_index ~= round((1-MAX_ZETA)*k)) && (up_index ~= 0 && up_index ~= min(round((1+MAX_ZETA)*k),length_Ykk-1))
                    if f_buf(low_index) ~= f_buf(k) && f_buf(up_index) ~= f_buf(k)
                        disp(['kBW: ' num2str(kBW) ', f_res-: ' num2str(f_buf(low_index)) ', f_res: ' num2str(f_buf(k)) ' Hz, f_res+: ' num2str(f_buf(up_index))]);
                        break;
                    end
                else
                    kBW = kBW - kBW_step;
                end
            end
            %% considering very small curve-drop in one side 
            if kBW <= 0 && (low_index_kBWmax ~= 0 || up_index_kBWmax ~= 0) && abs(Zm_buf(k,m)) > 5* ZM_LOW_LIMIT
                disp(['low_index_kBWmax: ' num2str(low_index_kBWmax) ', up_index_kBWmax: ' num2str(up_index_kBWmax)]);
                kBW = kBW_max;
                if low_index_kBWmax == 0
                    low_index = k;
                else
                    low_index = low_index_kBWmax;
                end
                up_index = max(k,up_index_kBWmax);
                disp(['kBW: ' num2str(kBW) ', f_res-: ' num2str(f_buf(low_index)) ', f_res: ' num2str(f_buf(k)) ' Hz, f_res+: ' num2str(f_buf(up_index))]);
            end
            if kBW <= 0
                disp(['Mode at ' num2str(f_buf(k)) ' Hz without both sides xAk !!!------------------------------']);
                %disp(['Zm_no=' num2str(m) ', |Zm|=' num2str(Zm_buf(k,m))]);
                continue;
            end

            lower_index = low_index;
            upper_index = up_index;
            f_upper= f_buf(upper_index);
            f_lower = f_buf(lower_index);
            %disp(['f_buf(lower_index): ' num2str(f_buf(lower_index)) ', f_buf(upper_index): ' num2str(f_buf(upper_index))]);
            % if (lower_index == 2 && upper_index == up_index) || (f_upper == f_buf(k) && f_lower == f_buf(k))
            %     disp(['Step 2: Please check mode: ' num2str(f_buf(k)) ' Hz------------------------------']);
            %     continue;
            % end
               
            if 1
                %% compensation: quite effective for zeta from 0.2 to 0.4, however not so for more higher zeta
                if f_upper == f_buf(k) && f_lower ~= f_buf(k)
                    bandwidth = 2*(f_buf(k)-f_lower);
                    BW_index=1;
                elseif f_lower == f_buf(k) && f_upper ~= f_buf(k)
                    bandwidth = 2*(f_upper-f_buf(k));
                    BW_index=3;
                elseif abs(f_upper+f_lower-2*f_buf(k))/f_buf(k) > 0.02
                    %% when asymmetry level is above a limit
                    BW_vec = [2*(f_buf(k)-f_lower)     f_upper-f_lower    2*(f_upper-f_buf(k))];
                    [bandwidth BW_index] = min(BW_vec(BW_vec>0));
                else
                    bandwidth = f_upper-f_lower;
                    BW_index=2;
                end
            end
            %disp(['wk- = ' num2str(f_lower) ' Hz, wk = ' num2str(f_buf(k)) ' Hz, wk+ = ' num2str(f_upper) ' Hz']);
            disp(['Bandwidth variant: ' num2str(BW_index)]);
            % Calculate Q factor
            %Q_factor = f_buf(k) / bandwidth;
            %zeta_mag = 1 / 2 / Q_factor;
            zeta_mag = bandwidth / (2*sqrt(1/kBW^2-1)*f_buf(k));

            %% Optimization based on std-2nd-order-system (not standard series RLC resonance circuit )
            % if BW_index == 1
            %     zeta_mag = sqrt( 0.5 - 0.5./sqrt((1-(1-zeta_mag).^2).^2+1) );
            % elseif BW_index == 3
            %     zeta_mag = sqrt( 0.5 - 0.5./sqrt(((zeta_mag+1).^2-1).^2+1) );
            % else
            %     zeta_mag = zeta_mag;
            % end
            %% The following compensation works for the IEEE 14 bus system
            % if any(zeta_mag <= DAMPING_UP_LIMIT) && any(zeta_mag > 0.1)
            %     zeta_mag = zeta_mag / sqrt(1-4*zeta_mag^2);
            % end

            %% Limiting the concerned damping range
            if zeta_mag >= DAMPING_UP_LIMIT %1/sqrt(2)
                disp(['Mode with large damping: ' num2str(zeta_mag) ', neglected!']);
                continue;
            end

            if 1
            zeta_sign0 = sign(real(Zm_buf(k,m))); %% Passivity
            %zeta_sign0 = -sign( real(Zm_buf(k,m)) * imag(Zm_buf(k,m)) ); %% PMD
            %% Still needs to distinguish zero-crossings for Rm or Xm and identify their relationship with damping sign
            dRm_dw_sign = sign(real(Zm_buf(k+1,m))-real(Zm_buf(k,m)));
            dXm_dw_sign = sign(imag(Zm_buf(k+1,m))-imag(Zm_buf(k,m)));

            if dRm_dw_sign == dXm_dw_sign %% possibly having two maximums
                if abs(real(Zm_buf(k,m))) > abs(real(Zm_buf(k-1,m))) && abs(real(Zm_buf(k,m))) > abs(real(Zm_buf(k+1,m))) %% at maximum
                    if real(Zm_buf(k,m)) * real(Zm_buf(k-1,m)) < 0
                        dRm_dw_sign = sign(real(Zm_buf(k,m))-real(Zm_buf(k-1,m)));
                    elseif real(Zm_buf(k+1,m)) * real(Zm_buf(k,m)) < 0
                        dRm_dw_sign = sign(real(Zm_buf(k+1,m))-real(Zm_buf(k,m)));
                    elseif abs(real(Zm_buf(k,m))-real(Zm_buf(k-1,m))) > abs(real(Zm_buf(k,m))-real(Zm_buf(k+1,m)))
                        dRm_dw_sign = sign(real(Zm_buf(k,m))-real(Zm_buf(k-1,m)));
                    else
                        dRm_dw_sign = sign(real(Zm_buf(k+1,m))-real(Zm_buf(k,m)));
                    end
                end
                if abs(imag(Zm_buf(k,m))) > abs(imag(Zm_buf(k-1,m))) && abs(imag(Zm_buf(k,m))) > abs(imag(Zm_buf(k+1,m))) %% at maximum
                    if imag(Zm_buf(k,m)) * imag(Zm_buf(k-1,m)) < 0
                        dXm_dw_sign = sign(imag(Zm_buf(k,m))-imag(Zm_buf(k-1,m)));
                    elseif imag(Zm_buf(k+1,m)) * imag(Zm_buf(k,m)) < 0
                        dXm_dw_sign = sign(imag(Zm_buf(k+1,m))-imag(Zm_buf(k,m)));
                    elseif abs(imag(Zm_buf(k,m))-imag(Zm_buf(k-1,m))) > abs(imag(Zm_buf(k,m))-imag(Zm_buf(k+1,m)))
                        dXm_dw_sign = sign(imag(Zm_buf(k,m))-imag(Zm_buf(k-1,m)));
                    else
                        dXm_dw_sign = sign(imag(Zm_buf(k+1,m))-imag(Zm_buf(k,m)));
                    end
                end
            end

            zeta_sign0x = -sign( real(Zm_buf(k,m)) ) * dXm_dw_sign; %% PMD
            zeta_sign0r = sign( imag(Zm_buf(k,m)) ) * dRm_dw_sign; %% PMD

            disp(['Rm = ' num2str(real(Zm_buf(k,m))) ', Xm = ' num2str(imag(Zm_buf(k,m)))] );
            disp(['sign of dXm/dw = ' num2str(dXm_dw_sign) ', sign of dRm/dw = ' num2str(dRm_dw_sign) ] );
            disp(['zeta_sign0x=-Rm*dXm/dw: ' num2str(zeta_sign0x) ', zeta_sign0r=Xm*dRm/dw: ' num2str(zeta_sign0r)]);

            %% mutual validation: if not equal, then throw the result
            if zeta_sign0x ~= zeta_sign0r ...
                    && dXm_dw_sign ~= 0 && dRm_dw_sign ~= 0 && real(Zm_buf(k,m))~=0 && imag(Zm_buf(k,m))~=0
                kRX=4; %% a threshod Rm/Xm ratio value for mutual validation
                if abs(real(Zm_buf(k,m))) >= abs(imag(Zm_buf(k,m)))/kRX && abs(real(Zm_buf(k,m))) <= kRX*abs(imag(Zm_buf(k,m)))
                    disp(['Mode at ' num2str(f_buf(k)) ' Hz skipped (damping sign mutual check failed) !!!------------------------------']);
                    continue;
                end
            end

            end

            %% Using impedance ratio (Zm_wn/Zm_peak) to determine the sign of the damping
            % %f_n = f_buf(k) * 1/sqrt(1-2*zeta_mag^2);
            f_n = f_buf(k) * 1/sqrt(sqrt(1+8*zeta_mag^2)-4*zeta_mag^2);
            [~, fn_index] = min(abs(f_buf(k:end))-f_n);

            %% Using impedance ratio (Zm_pole/Zm_peak) to determine the sign of the damping
            f_pole = f_buf(k);
            %f_pole = f_buf(k) * sqrt(1-zeta_mag^2)/sqrt(1-2*zeta_mag^2);
            %% for standard RLC-Resonance-circuit
            f_pole = f_buf(k) * sqrt(1-zeta_mag^2)/sqrt(sqrt(1+8*zeta_mag.^2)-4.*zeta_mag.^2);

            % disp(['Pole frequency: ' num2str(f_pole)]);
            [~, fp_index] = min(abs(f_buf(k:end)-f_pole));

            if zeta_mag <= DAMPING_UP_LIMIT
                if abs(imag(Zm_buf(k,m))) < abs(real(Zm_buf(k,m))) 
                    zeta_sign = zeta_sign0x;
                else
                    zeta_sign = zeta_sign0r;
                end
                %zeta_sign = zeta_sign_dtheta;
            else%if zeta_mag <= 0.05
                zeta_sign = zeta_sign1;
                %zeta_sign = zeta_sign2;
            end
            Damping = zeta_sign * zeta_mag;
            %% using the modal impedance at the actual oscillation frequency, i.e.wp, or natural frequency wn to calculate the sign of the damping
            
            disp(['Zm peak Frequency: ' num2str(f_buf(k)) ' Hz']);
            %disp(['----------- Est. damped Frequency: ' num2str(f_buf(fp_index+k-1)) ' Hz-----------']);
            %disp(['Bandwidth: ' num2str(bandwidth) ' Hz']);
            %disp(['Q Factor: ' num2str(Q_factor)]);
            %disp(['Est. natural Frequency: ' num2str(f_buf(fn_index+k-1)) ' Hz']);
            disp(['---------------Damping: ' num2str(Damping) '--------------']);
            
            res_no=res_no+1;
            ResZm(res_no)=Zm_buf(k,m);
            ResYm(res_no)=Ym_buf(k,m);
            ResFm(res_no) = f_buf(k);%k*Hstep;
            Mode_buf(res_no,:)=[k Zm_buf(k,m) Ym_buf(k,m) m f_buf(k) Damping f_buf(fp_index+k-1)];
        end
    end
end
%ResFm,abs(ResZm),
if res_no == 0
    display('No oscillation risk identified!!!');
    RMA_out = zeros(2,6);
    return,
end
RMA_out = Mode_buf;

if IF_PLOT_PEAK_FLAG== 1
%% plotting resonance peaks
subplot(3,1,1)
for k=1:length(ResFm)
    if abs(ResZm(k))< ZM_LOW_LIMIT
        continue;
    end
    %pp=plot(ResFm(k)*50,abs(ResZm(k)),'k.');
    pp=plot(ResFm(k),abs(ResZm(k)),'k.');
    pp.MarkerSize = 9;
    if ResFm(k) < 50 && ResFm(k) > 0.1
        txt = sprintf('%4.3fHz',ResFm(k));
    else
        txt = sprintf('%4.0fHz',ResFm(k));
    end
    %tt=text(ResFm(k)*50,abs(ResZm(k))+1,txt);
    tt=text(ResFm(k),abs(ResZm(k)),txt);
    tt.FontSize=9;
    tt.Color='r';
end
% subplot(3,1,2)
% for k=1:length(ResFm)
%     if abs(ResZm(k))< ZM_LOW_LIMIT
%         continue;
%     end
%     %pp=plot(ResFm(k)*50,abs(ResZm(k)),'k.');
%     pp=plot(ResFm(k),real(ResZm(k)),'k.');
%     pp.MarkerSize = 9;
%     if ResFm(k) < 50 && ResFm(k) > 0.1
%         txt = sprintf('%4.3fHz',ResFm(k));
%     else
%         txt = sprintf('%4.0fHz',ResFm(k));
%     end
%     %tt=text(ResFm(k)*50,abs(ResZm(k))+1,txt);
%     tt=text(ResFm(k),real(ResZm(k)),txt);
%     tt.FontSize=9;
%     tt.Color='r';
% end
subplot(3,1,3)
for k=1:length(ResFm)
    if abs(ResZm(k))< ZM_LOW_LIMIT
        continue;
    end
    %pp=plot(ResFm(k)*50,abs(ResZm(k)),'k.');
    pp=plot(ResFm(k),imag(ResZm(k)),'k.');
    pp.MarkerSize = 9;
    if ResFm(k) < 50 && ResFm(k) > 0.1
        txt = sprintf('%4.3fHz',ResFm(k));
    else
        txt = sprintf('%4.0fHz',ResFm(k));
    end
    %tt=text(ResFm(k)*50,abs(ResZm(k))+1,txt);
    tt=text(ResFm(k),imag(ResZm(k)),txt);
    tt.FontSize=9;
    tt.Color='r';
end
end
%% Participation Factor Analysis
for k=1:length(ResFm)
    f=ResFm(k);
    Ykk_index=find(f_buf==f,1);% or use floor()
    Ykk_res_array=Ykk_buf(Ykk_index,:); 
    Ykk_res = reshape(Ykk_res_array,[nbus,nbus]).';
    Zkk=inv(Ykk_res);
    [T,Ym,L]=eig(Ykk_res);
    Zm_diag=1./diag(Ym);
    %% MODIFIED on 20.11.2023
    Zdiff=abs(Zm_diag-ResZm(k));
    [Min_diff, grp_no] = min(Zdiff);
    %% END OF MODIFICATION
    T=inv(L);
    PFkk=L(:,grp_no)*T(grp_no,:);%PF(N,N)represents participation factor of bus N 
    PFk=abs(diag(PFkk)).';
    [Max_Zm, Index_m] = max(abs(ResZm));
    if abs(ResZm(k)) == Max_Zm
        %fprintf('fm=%.2fHz\nZm = %.2f * exp(%.2f°j)\n',f,abs(ResZm(k)),rad2deg(angle(ResZm(k)))),
        %fprintf('Ym = %.6f * exp(%.2f°j)\n',abs(ResYm(k)),rad2deg(angle(ResYm(k)))),
        CRITICAL_MODE_NO = Mode_buf(Index_m,4);
        %format bank %shortG
        %format compact
        %abs(PF),
    end
    %Ykk_index, grp_no,
    %PFtable(:,k)=[ResFm(k)*50 abs(ResZm(k)) abs(diag(PF)).'];%PFkk(1,1) PFkk(1,1) PFkk(1,1)];%%Participation factor table
    %format bank
    PFtable(:,k)=[ResFm(k) abs(ResZm(k)) Mode_buf(k,6) PFk];
    %PFtable(:,k)=[ResFm(k) ResZm(k) abs(diag(PF)).'];%PFkk(1,1) PFkk(1,1) PFkk(1,1)];%%Participation factor table
end

%% GET PF CURVES for the CRITICAL MODE obtained from the Ym or the Zm envelope curve
for k=1:1:length(f_buf)
    Ykk_m_array=Ykk_buf(k,:); 
    Ykk_mk = reshape(Ykk_m_array,[nbus,nbus]).';
    Zkk=inv(Ykk_mk);
    [T,Ym,L]=eig(Ykk_mk);
    %% the diag(Ym) is not necessarily equal to Ym_buf(k,:).' since they can be switched at intersecting points
    %% thus the PFs to the critical modes at critical frequency areas still need to be checked
    %abs(1./diag(Ym))
    [Min_Ym, grp_no] = min(abs(diag(Ym)));
    Ym_diag_diff = Min_Ym - Ym_buf(k,:);
    [Min_diff, MODE_NO] = min(Ym_diag_diff);
    %fprintf('|Zm%0.0f|=%4.2f, grp_no==%0.0f, f=%4.2fHz\n',MODE_NO,1/Min_Ym,grp_no,f_buf(k));
    T=inv(L);
    PFkk=L(:,grp_no)*T(grp_no,:);%PF(N,N)represents participation factor of bus N
    PFk=abs(diag(PFkk)).'; % PFk is a vector while PFkk is a matrix
    PF_buf(k,:) = PFk;
end

%% Plotting participation curves
if IF_PLOT_PF_FLAG == 1
    figure,
    PORT_NUM = nbus;
    if PORT_NUM == 1
    semilogx(f_buf,PF_buf(:,1)),
    legend('dc');
    elseif PORT_NUM == 2
    semilogx(f_buf,PF_buf(:,1), f_buf, PF_buf(:,2)),
    legend('p','n')
    elseif PORT_NUM == 3
    semilogx(f_buf,PF_buf(:,1),line_style(1), f_buf, PF_buf(:,2),line_style(2),f_buf, PF_buf(:,3),line_style(3)),
    legend('\itp\rm_P','\itp\rm_N','\itp\rm_D_C')
    elseif PORT_NUM == 4
    semilogx(f_buf, PF_buf(:,1), f_buf, PF_buf(:,2), f_buf, PF_buf(:,3), f_buf, PF_buf(:,4)),
    legend('N1','N2','N3','N4')
    elseif PORT_NUM == 6
    semilogx(f_buf, PF_buf(:,1),'r', f_buf, PF_buf(:,2),'m',f_buf, PF_buf(:,3),'b', f_buf, PF_buf(:,4),'--b',f_buf, PF_buf(:,5),'--r',f_buf, PF_buf(:,6),'--m'),
    legend('MMC1-P','MMC1-N','MMC1-DC','MMC2-DC','MMC2-P','MMC2-N')
    elseif PORT_NUM == 7
    semilogx(f_buf, PF_buf(:,1), f_buf, PF_buf(:,2),f_buf, PF_buf(:,3), f_buf, PF_buf(:,4),f_buf, PF_buf(:,5),f_buf, PF_buf(:,6),f_buf, PF_buf(:,7)),
    legend('N1','N2','N3','N4','N5','N6','N7')
    end
    hold on,grid on,
    legend boxoff;
    title('Participation Factors (PFs)');% in Ohms'),
    ylabel('Magnitude');% in Ohms'),
    xlabel('{\it f} in Hz'),
    % xticks([0.1 1 10 100 1000]),xlim([1 200]),
    % xticklabels({'0.1','1','10','100','1000'})
end

%fprintf('Participation factors at frequencies of modal resonance peaks:\n');
%format bank %shortG
%ResZm=abs(ResZm),
PFtable; %Table should also give resonance mode no
%END of PF calculation

%%%%%%%%%%%%%%%%%%%%%CALCULATION OF SENSITIVITY INDEX%%%%%%%%%%%
%%%%NOT INCLUDED YET, REFERENCE SCRIPT IS <test3bus_Ref for RMA PF and SI>