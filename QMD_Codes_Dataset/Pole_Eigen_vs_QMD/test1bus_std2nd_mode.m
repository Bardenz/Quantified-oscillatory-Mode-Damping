%% symbolic analysis

if 0
syms w s D_sym wn_sym wn zeta Ro(w) Xo(w) D GoM(w) R L
Ykk_sym = ((j*w)^2 + 2*D_sym*wn_sym*j*w + wn_sym^2) / wn_sym^2;
Zm_sym = inv(Ykk_sym);
Rm_sym = wn_sym^2*(wn_sym^2-w^2) / ((wn_sym^2-w^2)^2 + (2*D_sym*wn_sym*w)^2);
Xm_sym = -2*D_sym*wn_sym^3*w / ((wn_sym^2-w^2)^2 + (2*D_sym*wn_sym*w)^2);
Xm_sym = (wn_sym^2-w^2)^2 + (2*D_sym*wn_sym*w)^2;
dX_dw = diff(Xm_sym,w);

Rmw = ( (wn^2-w^2)*Ro(w)+2*D*wn*w*Xo(w) )*wn^2 / ((wn^2-w^2)^2 + (2*D*wn*w)^2);
Xmw = (-2*D*wn*w*Ro(w)+(wn^2-w^2)*Xo(w) )*wn^2 / ((wn^2-w^2)^2 + (2*D*wn*w)^2);

dRmw_dw = diff(Rmw,w),
dXmw_dw = diff(Xmw,w),

YmMag = ((wn^2-w^2)^2+(2*D*wn*w)^2)/wn^2/GoM(w);
ZmMag = sqrt(R^2+w^2*L^2) *wn^2 / sqrt((wn^2-w^2)^2 + (2*D*wn*w)^2);
ZmMag = sqrt(Ro(w)^2+Xo(w)^2) *wn^2 / sqrt((wn^2-w^2)^2 + (2*D*wn*w)^2);
%ZmMag = R*wn^2 / sqrt((wn^2-w^2)^2 + (2*D*wn*w)^2);
dYmMag_dw = diff(YmMag,w),
dZmMag_dw = diff(ZmMag,w),

poly = w^4 - (1-2*zeta^2)*wn^2*w^2 + (1-8*(1-zeta^2)*zeta^2)*wn^4;
solve(poly,w);

zeta = [0.1 0.2 0.3 0.4 0.5];
zeta_sim = [0.1 0.194 0.3 0.44 0.45];
wk_sim = [0.99 0.96 0.905 0.825 0.705];
wk_sim_zeta_ref = wk_sim .* sqrt(1-zeta.^2)./sqrt(1-2*zeta.^2);
wk_sim_zeta_sim = wk_sim .* sqrt(1-zeta_sim.^2)./sqrt(1-2*zeta_sim.^2);

return,

end

x = 1/sqrt(2);
%x = 0.99999;
c_zm = 4/(x^2);
zeta = -0.707:0.0001:0.707;
wk = sqrt( 1-2.*zeta.^2);
w_k1 = sqrt( 1-2.*zeta.^2 + abs(zeta).*sqrt((c_zm-4).*(1-zeta.^2)) );
%zeta2 = -0.3827:0.0001:0.3827;
w_k2 = sqrt( 1-2.*zeta.^2 - abs(zeta).*sqrt((c_zm-4).*(1-zeta.^2)) );
%size(zeta), size(w_k1), size(w_k2),

%subplot(1,2,1)
plot(zeta, w_k1, zeta, w_k2, zeta, 2*wk-w_k1,'--'),legend on
xlabel('{\it\zeta} in pu'), ylabel('\omega / \omega_n');
legend('\omega_k+','\omega_k-','\omega_k- (fitted)')

w_k22 = sqrt( 1-2.*zeta.^2 - 2.*abs(zeta).*sqrt(1-zeta.^2));
w_k22 = w_k2;
bw = (w_k1-w_k22);
bw_top = 2.*(w_k1-wk);
bw_btm = 2.*(wk-w_k22);
x=0;
for zeta_cur = -0.707:0.0001:0.707
    x=x+1;
    bw_cmp(x) = min([2*(w_k1(x)-wk(x)) 2*(wk(x)-w_k22(x)) (w_k1(x)-w_k22(x))]);
end
%size(zeta), size(bw), size(bw_cmp),

%subplot(1,2,2)
%% sign of zeta to consider
zeta_bw = sign(zeta).*bw./(2.*wk);
zeta_top =sign(zeta).*bw_top./(2.*wk);
zeta_btm = sign(zeta).*bw_btm./(2.*wk);
%zeta_cmp = sign(zeta).*bw_cmp./(2.*wk);
zeta_cmp = bw_cmp./(2.*wk);
zeta_cmp2 = sign(zeta).*sqrt( 0.5 - 0.5./sqrt(((zeta_cmp+1).^2-1).^2+1) );
%semilogy(zeta, zeta_zm1, zeta, zeta_zm2)
plot( zeta, zeta_bw, zeta, zeta_top, zeta, zeta_btm, zeta, zeta_cmp2), grid on
xlabel('{\it\zeta} in pu'), ylabel('{\it\zeta}_e_s_t in pu');
legend('\zeta_t_o_p_-_b_t_m','\zeta_i_n_i_t','\zeta_b_t_m','\zeta_b_w'), ylim([-0.65 0.65]), 

% subplot(1,3,3)
% [fittedcurve,gof] = fit(zeta.', zeta_zm2.', 'poly2')
% plot(fittedcurve, zeta.', zeta_zm2.'), ylim([-1 1]),
% legend('Location','NorthWest');


return,
zeta_vec = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
for k=1:1:length(zeta_vec)
    zeta = zeta_vec(k);
    fn = 0.01:0.01:1000;
    wn = 2*pi.*fn; % rad/s
    R_zeta = 1 ./ (4 - 3*zeta.^2);
    dXm_dw = -2./((4-3*zeta.^2).*zeta.*wn) + 8 * sqrt(1-zeta.^2) ./ ((4-3*zeta.^2).^2.*zeta.*wn.^2);
    zeta_sign = sign(-R_zeta .* dXm_dw) .* zeta;
    %plot(zeta,zeta_sign),grid on
    semilogx(fn,zeta_sign),grid on,hold on
end
xlabel('Natural frequency {\it f}_n in Hz'),
%xlabel('{\it\zeta} - Time-domain'),
ylabel(['Extracted {\it\zeta}']);
title('Considering only the sign of the extracted damping')
ylim([-2 2]),
legend('\zeta = 0.1','\zeta = 0.2','\zeta = 0.3','\zeta = 0.4','\zeta = 0.5','\zeta = 0.6','\zeta = 0.7','\zeta = 0.8','\zeta = 0.9')