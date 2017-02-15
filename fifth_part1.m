
%%  Program for Optimal Policy for RF energy Harvesting Rate
%%% Network model: Time Switching Architecture
%   That is, either the information receiver or the RF energy harvester is connected
%   to the antenna at a given time.



%%  Variable description
% N       : total number of transmitters
% K       : Point Process set
% tau     : portion of time to harvest energy
% beta    : RF to DC power conversion efficiency
% Ps      : transmit power of RF transmitter k
% gamma   : path loss exponent
% h       : channel power gain
% theta   : rate parameter
% delta   : shape parameter
% d       : distance between transmit antenna of an RF transmitter to the
%           reciever antenna of RF powered device
% x       : coordinates of k-th Rf transmitter referential to rf Device
% mu      : a constant
% max     : maximum value for Ps
% xie     : density of ambient RF transmitter
% apsi    : closest distance of transmitter to the device
% R       : Range of harvestor
% gamma_A : path loss exponent for access point
% Pa      : Power transmitted by access point
% da      : Distance of access point from harvester
% W       : Transmission Bandwidth (The b.w. of channel between the access point and the device)
% Pc      : A base circuit power device will consume(minimum power required)
% ic      : Interference Coefficient
% alpha   : Repulsion factor
% sigs    : Power density of AWGN  
% P_H     : RF energy harvesting rate(in watts) by the device from Rf
%           transmitter k
% P_HT    : Aggregated RF energy harvesting rate by the device equipped with
%           time switching architecture

% F,F_bar : Used to define thermal noise and noise in channel
% C_TS    : Maximum Transmission Rate
% Peo     : Energy outage probability
% K_xy    : Ginibre Kernel
% A       : Integral operator with kernel
%% Initialization and Input parameter settings
N       = 5;
tau     = 0:0.001:1;
lamb    = 25; %%% *%%%%%%%"check"*
delta   = 1; %so that 'hk' can be equal to 'exp(theta)'
beta    = 0.3; % conversion efficiency : 30%
mu      = 0.01; %must be grater than 0
e       = 0.05;
gamma   = 4;
gamma_A = 5;
theta_A = 20;
R       = 30;  % 30 m
Pa      = 39.81; %39.81W =46dBm
Ps      = 10; %100mW
lemda   = 2;
theta   = 0.5; %m/omega
da      = 1;  %(2.5 optimized)w
W       = 20*(10^6);
Pc      = 2.64*(10^(-6));
ic      = 0.5; %(ic->[0,1])
alpha   = 0; %aplha=0 for poisson
sigs    = 10^(-15); %sig^2 = -120dbm/hz = 1*10(-15) W/hz
xie     = [0.005 0.01];%[0 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1]; %needs to be variable
h       = zeros(N,1);
d       = zeros(N,1);
%% Input data 
% All data that needed to be entered by the user
%input('Portion of time to harvest energy: ');
%input('RF to DC power conversion efficiency');
%randi([0,MAX],[K,1]); %power value changes from '0' to 'max' with 'K'(no of transmitter)by 'T' (instances of time)
%input('Rate parameter: '); %theta must be positive, for all transmitters(1 by k)
% xk =[50 40 70; 60 80 90; 59 92 63; 50 94 133;120 50 31]; %input('Enter coordinates: '); %put array of coordinates for K no of transmitters,'(K,3)'
P_HTS=zeros(1,1001);
C_TS=zeros(1,1001);
si =zeros(2,1);
Peo=zeros(1,1001);
% Lx =zeros(N,1001);
% Ly =zeros(N,1001);
% G_k=zeros(N,2);
%% Defining different functions to be used further
%%% Distance of transmitter from device
for w=1:6 
    
  for j=1:1001
      x=R*random('poisson',lamb,[N,2]);% x is matrix with two columns for x and y and N rows for No. of transmitters
        
      for i=1:N
        d(i,1) = e+(((x(i,1))^2+x(i,2)^2)^0.5);   
      end
   
%%% Channel gain between transmitter and device
    
      h(:,1) = random('gamma',delta,theta,[N,1]);
    
%%% Channel gain between Access point and device
     ha  = lemda*exp(lemda*da);
     acc = (Pa*ha)/((da)^gamma_A);
%%% Channel and internal Noise approx. function
     F = mu*exp(-mu);
     u = mu*( -mu + 1 ); 
     syms v
     Ei = double(-( int((1/(mu*(v+1)))*exp(-mu*(v+1)),v,mu,inf, 'PrincipalValue',true))); % Calculate value of the integral in the paper mailed page no 9 and conver the result into double 
     F_bar = (1/(-mu*exp(mu)*Ei));

%% Defining Ginibre process
%     for i=1:N
%      G_K(i,k)   = xie(k)*(exp(pi*xie(k)*x(i,1)*x(i,2))*exp(-(pi/2)*(x(i,1)^2+x(i,2)^2))); %K(x,y)[number of row changes with transmitter, number of coulm changes with xie]
%     end
%     for j=1:1001
%         for i=1:N
%          Lx(i,j)  = ((1-((1+(mu*tau(j)*beta*Ps*F_bar)/(theta*Pc*((x(i,1)+e)^gamma)))^(-1)))^0.5);
%          Ly(i,j)  = ((1-((1+(mu*tau(j)*beta*Ps*F_bar)/(theta*Pc*((x(i,2)+e)^gamma)))^(-1)))^0.5);% number of rows changes with transmitter and number of columns changes with tau
%         end
%     end
%% Calculating Energy Harvesting rate
%%% RF energy harvesting rate instantaneous
     P_H = zeros(N,1);
    
      for i=1:N
      P_H(i) = (Ps*h(i,1))/((d(i))^gamma);
      end
   

%%%Interference of ambient tansmitter
     I_TS= sum(P_H);
%%% Aggregate RF energy harvesting rate
     s=I_TS+acc;
     
      P_HTS(1,j) = (tau(1,j)*beta)*((F_bar)/(1+F))*s;  % for the expression before 's' we can have some other expression to show noise


%%% Expectation of RF energy harvesting rate 
     cons = (Pa/((theta_A)*(da^gamma_A)));
     
     for k=1:2
      si(k)=(2*pi*(xie(k))*Ps*delta)/theta;
     end
     
     if (gamma~=1 && gamma~=2 )
       Ir = ((e^(2-gamma))-(((R+e)^(1-gamma))*(e+(gamma-1)*R)))/((gamma-2)*(gamma-1));
     else if (gamma==1)
       Ir = R -(e*log(1+(R/e)));
     else Ir = (log(1+(R/e)))- (R/(R+e));
         end
     end
     
     E_PHTS=zeros(2,1);
   
     for i=1:2
      E_PHTS(i)= (tau(1,j)*beta)*(cons+((si(i))*Ir));
     end
    
%% Performance Matrics     
%%% MAXIMUM TANSMISSION RATE OF ACCESS POINT
   
        if P_HTS(1,j) >= Pc    
         C_TS(1,j)=(1-tau(1,j))*W*(log(1+((acc)/((sigs)+(ic*I_TS)))));
        else
         C_TS=0;
        
    end
%% ENERGY OUTAGE PROBABILITY
    
       Peo(1,j)=double(((1+(mu*tau(1,j)*Pa*F_bar)/((theta_A)*Pc*(da^gamma_A)))^(-1)));
  end
   
%% PLOTS AND RESULTS
%%% Energy outage probability versus tau
%figure
%plot(tau,Peo);
% xlabel('Time switching coefficient, tau')
% ylabel('Energy outage probability')
% title('Energy outage probability versus Time switching coefficient,tau')

% %%% Performance Analysis versus tau

    [AX,H1,H2]=plotyy(tau,P_HTS,tau,C_TS,'semilogy','semilogy');
     
     set(get(AX(1),'Ylabel'),'String','Expectation of Energy Harvesting Rate')
     set(get(AX(2),'Ylabel'),'String','Maximum Transmission Rate')
     
%      xlabel('tau(portion of time for harvesting energy)')
%      title('Performance Analysis with tau') 
hold on
da=da+0.5;
end
%% Function for Fredholm determinant
% function d = DetNystrom(K,z,a,b,m)
% [w,x] = QuadratureRule(a,b,m);
% w = sqrt(w);
% [xj,xi] = meshgrid(x,x);
% d = det(eye(m)+z*(w'*w).*K(xi,xj));
% end
