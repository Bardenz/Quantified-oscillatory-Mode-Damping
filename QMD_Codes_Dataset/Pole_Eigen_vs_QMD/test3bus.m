% this is a test from the paper Harmonic Resonance Mode Analysis_Wilsun Xu
clear; close all; clc;
f1=50;%Hz
w1=2*pi*f1;
%nbus=3;
Hmin=0.02;
Hmax=20;
Hstep=1e-3;
%%%%TEST REPORT on 28.06.2018
%%%accuracy of calculation results depends on frequency step
%%%results approximately match with those of the Master Thesis:
%%%Harmonic Resonance Mode Analysis and Application for Offshore Wind Power Plants_2015
Td=1e-5;

for k=Hmin/Hstep:1:Hmax/Hstep
    h=k*Hstep;
    w=w1*h;

    B1=0.13j*h;
    Xsys=8+9j*h;
    X12=0.835+80j*h;
    X23=0.835+80j*h;
    B3=0.26j*h;

    % B1=0.0013j*h;
    % Xsys=0.04+0.3j*h;
    % X12=0.835+4j*h;
    % X23=0.835+4j*h;
    % B3=0.0013j*h;

    Gd = exp(-1.5*j*w*Td);
    % 
    % X23 = X23*Td;
    % X12 = X12*Td;
    %Xsys = Xsys*Td;
if 1       
    Ykk=[B1+1/Xsys+1/X12     -1/X12              0;
           -1/X12         1/X12+1/X23           -1/X23;
              0               -1/X23          B3+1/X23];
else %% lumped impedance along the path from Bus1 to bus3
    Z_INM=Xsys+1/(B1+1/(X12+X23+1/B3));
    Ykk=1/Z_INM;
end
    % [T,Ym,L]=eig(Ykk);
    % Zm_diag=(1./diag(Ym)).';
    % Zm_buf(k,:)=Zm_diag;
    % Zkk=inv(Ykk);  
    % Zkk_buf(k,:)=diag(Zkk).';
    Ykk_buf(k,:)=reshape(Ykk.',1,[]);
    f_buf(k)=h*f1;
end
RMAout = my_RMA_fun(f_buf,Ykk_buf,0);

%% symbolic analysis
syms s

B1_s=0.13*s/w1;
Xsys_s=8+9*s/w1;
X12_s=0.835+80*s/w1;
X23_s=0.835+80*s/w1;
B3_s=0.26*s/w1;

Ykk_sym=[B1_s+1/Xsys_s+1/X12_s     -1/X12_s              0
            -1/X12_s         1/X12_s+1/X23_s           -1/X23_s
                0               -1/X23_s          B3_s+1/X23_s];
Ykk_det = det(Ykk_sym);
Ykk_zeros = single( vpasolve(Ykk_det) );
fr=imag(Ykk_zeros)/pi/2,
damp_ratio = -real(Ykk_zeros)./abs(Ykk_zeros),
return,