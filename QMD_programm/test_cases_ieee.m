% Form Admittance matrix for IEEE-14 bus sysSTART_BUStems and perform
% Resonance Mode Analysis(RMA)
clear all

%ieee14bus_data
%case57_B60Hz
case118
%case14

%case118_eig_data

RMA_PLOT_FLAG = 0;

Hmin=2; Hmax=240; Hstep = 1/50; %% default for the tests in publication Hmin=10/50; Hstep = Hmin/100;

if START_BUS == 0
    BUS_OFFSET = 1;
else
    BUS_OFFSET = 0;
end

Cnt = 0;           
for k=Hmin/Hstep:1:Hmax/Hstep       % Frequency range
    h=k*Hstep;               % Harmonic order

    if h >= 1/50 && rem(h,0.5/50)~=0
        continue;
    end
    % if w >= 2*pi*1 && w < 2*pi*10 && rem(w,2*pi*0.1)~=0
    %     continue;
    % elseif w >= 2*pi*10 && w < 2*pi*100 && (rem(w,2*pi*1)~=0 || w==wb) %% Remove f1 component to prevent ERROR IN Zm frequency sweep
    %     continue;
    % elseif w >= 2*pi*100 && (rem(w,2*pi*5)~=0 || w==2*wb || w==3*wb) %% Remove f1 component to prevent ERROR IN Zm frequency sweep
    %     continue;
    % end
    Cnt = Cnt + 1;

    fb = branch_data(:,1) + BUS_OFFSET;  % From bus number...
    tb = branch_data(:,2) + BUS_OFFSET;  % To bus number...
    r = branch_data(:,3);   % Resistance, R...
    x = branch_data(:,4);   % Reactance, X...
    gb = branch_data(:,5);  % Ground Admittance,B
    z = r + j*x*h;          % Z matrix...
    y = 1./z;               % To get inverse of each element...
    b = j*gb*h;              % Make B imaginary... 

    Y = zeros(nbus,nbus);   % Initialise YBus...

    nbranch = length(fb);   % no. of branches...  
    % Formation of the Off Diagonal Elements...
    for xy=1:nbranch
        if BUS_OFFSET == 0
            if (tb(xy)==0) || (tb(xy)==301) || (tb(xy)==302)
                continue;
            end
        end
        %fprintf('off diagonal elements: from bus %d to bus %d\n',fb(k),tb(k));
        Y(fb(xy),tb(xy)) = Y(fb(xy),tb(xy))-y(xy);
        Y(tb(xy),fb(xy)) = Y(fb(xy),tb(xy));
    end
 
    % Formation of Diagonal Elements....
    for m =1:nbus
        for n =1:nbranch
            if BUS_OFFSET == 0
                if (m==9) && (fb(n)==9) && (tb(n)==0)%capacitor branch of bus 9      
                    Y(m,m) = Y(m,m) + b(n);
                elseif (m==8) && (fb(n)==8) && (tb(n)==0)%filter branches of bus 8
                    Y(m,m) = Y(m,m) + 1/(z(n)+1/b(n));
                elseif (m==3) && (fb(n)==3) && (tb(n)==0)%filter branches of bus 3
                    Y(m,m) = Y(m,m) + 1/(z(n)+1/b(n));
                elseif (m==3) && (fb(n)==3) && (tb(n)==301 || tb(n)==302)%converter transformer branches of bus 3
                    Y(m,m) = Y(m,m) + y(n);%continue;
                elseif (fb(n)==m) || (tb(n)==m)%series branches between buses    
                        % if b(n)~=0%long line correction,cable length is considered, but cancelled in the correcting factor
                        %     Z0=sqrt(z(n)/b(n));
                        %     r0=sqrt(z(n)*b(n));
                        %     %z(n)=z(n)*sinh(r0)/r0;
                        %     %b(n)=b(n)*tanh(r0)/(0.5*r0);
                        %     z(n)=Z0*sinh(r0);
                        %     b(n)=tanh(r0/2)/Z0;
                        %     y(n)=1/z(n);
                        %     b(n)=2*b(n);
                        % end
                    Y(m,m) = Y(m,m) + y(n) + b(n)/2;
                end
            else
                if (fb(n)==m) || (tb(n)==m)%series branches between buses    
                    Y(m,m) = Y(m,m) + y(n) + b(n)/2;
                end
            end
            %% End of n-loop
        end
        % end of m-loop
    end
 
    % Adding load impedance to admittance matrix
    ldN=bus_data(:,1) + BUS_OFFSET; % bus no
    if length(bus_data(1,:)) == 5
        ldP=bus_data(:,2); % 
        ldQ=bus_data(:,3);
        ldU=bus_data(:,4);
        ldUn=bus_data(:,5);
    elseif length(bus_data(1,:)) == 6
        ldP=bus_data(:,3); % 
        ldQ=bus_data(:,4);
        ldU=bus_data(:,5);
        ldUn=bus_data(:,6);
    end
    ldnum=length(ldN);
    for kk =1:ldnum
        if (ldP(kk)==0)%CONVENTIONAL LOAD MODEL
            Yld=1/(j*h*ldU(kk)^2/ldQ(kk));
        else%CIGRE LOAD MODEL: Type-C
            Rld=ldU(kk)^2/ldP(kk);
            %Yld=1/(Rld+j*h*0.073*Rld)+1/(j*h*Rld*(6.7*ldQ(kk)/ldP(kk)-0.74));
            %% update on 09.05.2025 to 
            Yld=1/(Rld+j*h*0.073*Rld)+1/(j*h*Rld/(6.7*ldQ(kk)/ldP(kk)-0.74));
        end
        no=ldN(kk);
        if BUS_OFFSET == 0
             if no==301 || no==302
                 no=3;
                 %Yld=1/(1/Yld+j*h*0.028);
             end
        end
        Y(no,no)=Y(no,no)+Yld;
    end
    
    % Adding generator impedance to admittance matrix
    for gen_no=1:1:length(gen_data(:,1))
        bus_no = gen_data(gen_no,1) + BUS_OFFSET;
        Y(bus_no,bus_no)=Y(bus_no,bus_no)+1/(0.0+gen_data(gen_no,6)*j*h);
    end
    
    %Y(3,3)=Y(3,3)+2/(j*0.01*h);
    
    % %CALCULATING MODAL IMPEDANCE
    % [T,Ym]=eig(Y);
    % zm=(1./diag(Ym))';
    % Zm(k,:)=zm;
    % 
    % Zkk=inv(Y);  
    % Zk(k,:)=diag(Zkk)';

    f_buf(Cnt)=h*50;
    Ykk_buf(Cnt,:)=reshape(Y.',1,[]);
end

RMAout = my_RMA_fun(f_buf,Ykk_buf,RMA_PLOT_FLAG);
if RMA_PLOT_FLAG == 1
    return,
end

%figure
fk_RMA = RMAout(:,5);
D_RMA = RMAout(:,6);
if 0
%fk_RMA2pole = RMAout(:,7);
%semilogx(fk_RMA, D_RMA,'rx'), hold on, grid on
plot(fk_RMA, D_RMA,'rx'), hold on, grid on
xlabel('{\it f} in Hz'), ylabel(['\zeta in pu']);
%title('Oscillatory modes')
legend('Eigenvalues (State-Space)','Modes (Z-Curves)')
%%legend('Poles (s-Domain)','Zeros (s-Domain)','Modes (Z-Curves)')
% xticks([1 10 100 1000 10000])
% xticklabels({'1','10','100','1000','10000'}),
else
sigma = -D_RMA.*fk_RMA*2*pi; wd = sqrt(1-D_RMA.^2).*fk_RMA*2*pi;
plot(sigma, wd,'rx'), hold on, grid on
legend('Eigenvalues (State-Space)','Eigenvalues (Z-Curves)')
end
