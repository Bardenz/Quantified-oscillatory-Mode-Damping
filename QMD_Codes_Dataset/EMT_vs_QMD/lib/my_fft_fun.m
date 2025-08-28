function fft_out = my_fft_fun(Meas_input,Ts_fft)

%my_MMC_param

Hmin=0.002;
Hmax=60;

F1=50;%Hz
T1=1/F1;
%Ts= Ts_Power; %% 20e-6;%
Fs=1/Ts_fft;
Num1s=1/Ts_fft; % data number of 1 second

% Tsimu at least 10 + Tstart (Tstart is the time for the system to be in steady state)

%%%%%TESTING REPORT%%%%%%%%%%
%%%Power spectral density method has helped improving the impedance phase
%%%curve, but not on the magnitude curve
%%%for the dfig_5mw_model_machine.mat file, set index_no=6, Tstart=16

    nn=1;

    %% ------------ Time Setting -----------------
    Tstop = 19;
    Twindow = 2;
    %% ----------------------------------------------
    L=Twindow/Ts_fft;

    Uabc1=Meas_input;
    %%%%VOLTAGE FFT TRANSFORM
    Ut1=Uabc1(Tstop*Num1s-L+1:Tstop*Num1s,:);
    Us1 = fft(Ut1);
    Us1= Us1(1:L/2+1,:)/L;
    Us1(2:end-1,:) = 2*Us1(2:end-1,:);

    %figure,
    f = Fs*(0:(L/2))/L;
    h=Fs/F1*(0:(L/2))/L;
    %h=Fs*(0:(L/2))/L;%%%subsynchronous
    f_end=length(h);
        
        Us1_buf=zeros(f_end,2);
        Us1_buf(:,1)=Us1(:,1);
if 0    
    figure,
    semilogx(f,abs(Us1(:,1))),%,f,abs(Us2_buf(:,1))*100,'-.'),
    hold on, grid on,
    xlabel('Frequency in Hz'),
    ylabel('U_h / U_1'),
    xlim([Hmin Hmax]*50),%ylim([1e-6 inf]),
    %set(gca,'xticklabel',[]),
    title('Voltage spectrum');
end
    
for f_no=1:length(f)
    if f(f_no) == 0
        fft_out.f0 = Us1(f_no,:).';
    elseif f(f_no) == F1
        fft_out.f1 = Us1(f_no,:).';
    elseif f(f_no) == 2*F1
        fft_out.f2 = Us1(f_no,:).';
    elseif f(f_no) == 3*F1
        fft_out.f3 = Us1(f_no,:).';
    end
end