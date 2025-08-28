%% this is a test from the paper
%% "Frequency-Domain Modal Analysis of the Oscillatory Stability of Power Systems With High-Penetration Renewables"
clear; close all; clc;
Hmin=0.01/50; Hmax=40; Hstep=Hmin/2;

fn_vec = [0.1 1 10 100 1000]; % nature frequency vector
fn_vec = [1000]; % nature frequency vector
%D_vec = [-0.3 -0.2 -0.1 -0.07 -0.05 -0.03 -0.01 0.01 0.03 0.05 0.07 0.1 0.2 0.3];
%D_vec = [-0.01 -0.05 -0.02 0.01 0.05 0.08 0.1];
D_vec = [0.1 -0.1];
%D_vec = [-0.1 -0.2 -0.3 -0.4 -0.5];
%D_vec = [0 0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];

%% Frequency-Loop
for fn_index =1:length(fn_vec)
    Hmin=fn_vec(fn_index)/50/10; Hmax=fn_vec(fn_index)/50*10; Hstep=Hmin/100;
        if fn_vec(fn_index) < 1
            Hmin=fn_vec(fn_index)/50/10; Hmax=fn_vec(fn_index)/50*10; Hstep=Hmin/1000;
        end
    for D_index = 1:length(D_vec) %% Damping-LOOP
        f_n = fn_vec(fn_index); % Hz
        w_n = 2*pi*f_n; % natural angular frequency
        D = D_vec(D_index); % Damping ratio/factor
        f_k = w_n * sqrt(1-D^2) / 2 /pi;
        %disp(['Loop fn = ' num2str(f_n) ', D = ' num2str(D)] );
        Ykk_buf=[]; f_buf=[];

for k=Hmin/Hstep:1:Hmax/Hstep
    h=k*Hstep; f1=50; w1=2*pi*f1; w=w1*h;
    %% MOST IMPORTANT VARIATIONS TO CONSIDER: zeta and Go_fun
    Ro = 20; Xo = 20/w1*w; %% Xo>>Ro: Zero-crossing of X at wn, Ro>>Xo: Zero-crossing of R at wn
    Go_fun = (Ro+j*Xo); %% ((j*w)^2/w_n^2+1); (2*D*w_n+1j*w); 
    Go_fun = (2*D*w_n+1j*w*(1+4*D^2));
    Ykk = ((j*w)^2 + 2*D*w_n*(j*w) + w_n^2) / w_n^2 / Go_fun;
    % R=0.0; %w_n=1;
    % Ykk = ((j*w)^2 + 2*D*w_n*(j*w) + w_n^2) / w_n^2 / ((j*w)^2*(2*D-R) + (j*w)*((2*D-R)*R+1) + R);

    Ykk_buf(k,:)=reshape(Ykk.',1,[]);
    f_buf(k)=h*f1;
end

RMA_PLOT_SIGN = 1;
RMAout = my_RMA_fun(f_buf,Ykk_buf, RMA_PLOT_SIGN);

if RMA_PLOT_SIGN == 1
    continue;
end

fk_RMA = RMAout(1,5);
D_RMA = RMAout(1,6);
%fk_RMA2pole = RMAout(1,7);
line_style=['c+' 'cx' 'co'
            'r+' 'rx' 'ro'
            'b+' 'bx' 'bo'
            'm+' 'mx' 'mo'
            'g+' 'gx' 'go'
            'k+' 'kx' 'ko'
            'r+' 'rx' 'ro'
            'b+' 'bx' 'bo'
            'm+' 'mx' 'mo'
            'g+' 'gx' 'go'
            'k+' 'kx' 'ko'
            'r+' 'rx' 'ro'
            'b+' 'bx' 'bo'
            'm+' 'mx' 'mo'
            'g+' 'gx' 'go'
            'k+' 'kx' 'ko'];
% +   natural frequency and damping from eigenvalues / poles
% x   damped frequency and damping from eigenvalues / poles
% o   oscillation frequency and damping from modal impedance peaks
if 1
%semilogx(f_n, D,'b+', f_k, D,'bx', fk_RMA, D_RMA, 'bo'), hold on, grid on
%semilogx(f_n, D,line_style(D_index,1:2), f_k, D,line_style(D_index,3:4), fk_RMA, D_RMA, line_style(D_index,5:6)'), hold on, grid on
semilogx(f_n, D,line_style(D_index,3:4), fk_RMA, D_RMA, line_style(D_index,5:6)'), hold on, grid on
xlabel('{\it f} in Hz'), ylabel(['Damping ratio in pu']);
% xticks([0.01 0.1 1 10 100 1000 10000])
% xticklabels({'0.01','0.1','1','10','100','1000','10000'}),
% xlim([0.1 10]),
% %ylim([-0. 0.3]),
end

continue;

%% symbolic analysis
syms w s D_sym wn_sym
Ykk_sym = (s^2 + 2*D*w_n*s + w_n^2) / w_n^2;
Ykk_sym_para = ((j*w)^2 + 2*D_sym*wn_sym*j*w + wn_sym^2) / wn_sym^2;
Rm_sym = wn_sym^2*(wn_sym^2-w^2) / ((wn_sym^2-w^2)^2 + (2*D_sym*wn_sym*w)^2);
Xm_sym = -2*D_sym*wn_sym^3*w / ((wn_sym^2-w^2)^2 + (2*D_sym*wn_sym*w)^2);
dX_dw = diff(Xm_sym,w),

%Ykk_sym = (s^2 + 2*D*w_n*s + w_n^2) / w_n^2 / (2*D*w_n+s);
Ykk_det = det(Ykk_sym);
Ykk_zeros = single( vpasolve(Ykk_det) );

fr=imag(Ykk_zeros)/pi/2,
damp_ratio = -real(Ykk_zeros)./abs(Ykk_zeros),
%% symbolic
%Zm_sym = inv (Ykk_sym);
% %% from symbolic to transfer functions constructed by coefficients 
% [nm,dn] = numden(Zm_sym);
% nmd = sym2poly(nm);
% dnd = sym2poly(dn);
% Zm_TF= tf(nmd,dnd);
% step(Zm_TF,3), hold on, grid on,
% legend('\zeta = 0','\zeta = 0.01','\zeta = 0.05','\zeta = 0.1','\zeta = 0.2','\zeta = 0.3','\zeta = 0.4','\zeta = 0.5','\zeta = 0.6','\zeta = 0.7','\zeta = 0.8','\zeta = 0.9','\zeta = 1'),
%legend boxoff;
%xlabel('{\it t} in s'),

% figure,
% bode(Zm_TF)

    end %% D loop - variation of damping ratio
end % fn loop - variation of oscillation frequency





