% Form Admittance matrix for IEEE-14 bus systems and perform
% Resonance Mode Analysis(RMA)
clear all
close all

%ieee14bus_data
case57
%case39
%case14

if START_BUS == 0
    BUS_OFFSET = 1;
else
    BUS_OFFSET = 0;
end
        
syms s Y Ykk_sym Ykk_det

    fb = branch_data(:,1) + BUS_OFFSET;  % From bus number...
    tb = branch_data(:,2) + BUS_OFFSET;  % To bus number...
    r = branch_data(:,3);   % Resistance, R...
    x = branch_data(:,4);   % Reactance, X...
    gb = branch_data(:,5);  % Ground Admittance,B
    %z = r + j*x*h;          % Z matrix...
    z = r + s*x/wb;         % Z matrix...
    y = 1./z;               % To get inverse of each element...
    b = s*gb/wb;             % Make B imaginary...

    nbranch = length(fb);   % no. of branches...
    Yval = zeros(nbus,nbus);   % Initialise YBus...
    %% Create a matrix of symbolic numbers
    Y = sym(Yval);

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
    [row_num,~] = size(bus_data);
    if row_num ~= 0
        ldN=bus_data(:,1)  + BUS_OFFSET; % bus no
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
            if ldP(kk)==0 %CONVENTIONAL LOAD MODEL
                Yld=1/(s/w0*ldU(kk)^2/ldQ(kk));
            else%CIGRE LOAD MODEL
                Rld=ldU(kk)^2/ldP(kk);
                %Yld=1/(Rld+s/wb*0.073*Rld)+1/(s/wb*Rld*(6.7*ldQ(kk)/ldP(kk)-0.74));
                %% update on 09.05.2025 to 
                Yld=1/(Rld+s/wb*0.073*Rld)+1/(s/wb*Rld/(6.7*ldQ(kk)/ldP(kk)-0.74));
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
    end

    % Adding generator impedance to admittance matrix
    [row_num,~] = size(gen_data);
    if row_num ~= 0
        for gen_no=1:1:length(gen_data(:,1))
            bus_no = gen_data(gen_no,1) + BUS_OFFSET;
            Y(bus_no,bus_no)=Y(bus_no,bus_no)+1/(0.0+gen_data(gen_no,6)*s/wb);
        end
    end
    
    %Y(3,3)=Y(3,3)+2/(j*0.01*h);
    
%% solving CLTF poles
if 1
    Ykk_sym = Y;
    Ykk_det = det(Ykk_sym);
else
    %Ykk_det = prod(diag(Y));
    lambda = eig(Y),
    [T,E] = eig(Y),
    Ykk_det = prod(lambda);
end

Ykk_zeros = single( vpasolve(Ykk_det) );
poles = Ykk_zeros(imag(Ykk_zeros)>0); %% get conjugate complex poles
fd=imag(poles)/pi/2;%% damped frequencies of poles
fn=abs(poles)/pi/2;
damping = -real(poles)./abs(poles); %% = -(-zeta*wn) / wn

if 1
    semilogx(fn,damping,'b+'), hold on, grid on
    %semilogx(fn,damping,'b+', fd,damping,'bsquare'), hold on, grid on
else
    %% solving CLTF zeros
    Zkk_zeros = single( vpasolve(inv(Ykk_det)) );
    zeros = Zkk_zeros(imag(Zkk_zeros)>0);
    fd_z=imag(zeros)/pi/2; %% damped frequencies of zeros
    damping_z = -real(zeros)./abs(zeros);
    semilogx(fn,damping,'b+',fd_z, damping_z,'bo'), hold on, grid on
end
%xlabel('{\it f} in Hz'), ylabel(['\zeta in pu']);
%test_cases_ieee