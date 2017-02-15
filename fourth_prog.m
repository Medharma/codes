
%%  Program for Optimal Policy for RF energy Harvesting Rate
%%% Network model: Time Switching Architecture
%   That is, either the information receiver or the RF energy harvester is connected
%   to the antenna at a given time.



%%  Variable description
% tau     : portion of time to harvest energy
% beta    : RF to DC power conversion efficiency
% Ps      : transmit power of RF transmitter k
% gamma   : path loss exponent
% hk      : channel power gain
% theta   : rate parameter
% delta   : shape parameter
% dk      : distance between transmit antenna of an RF transmitter to the
%           reciever antenna of RF powered device
% xk      : coordinates of k-th Rf transmitter referential to rf Device
% mu      : a constant
% max     : maximum value for Ps
% xie     : density of ambient RF transmitter
% apsi    : closest distance of transmitter to the device
% R       : Range of harvestor
% gamma_A : path loss exponent for access point
% Pa      : Power transmitted by access point
% da      : Distance of access point from harvester
% 
% P_H     : RF energy harvesting rate(in watts) by the device from Rf
%           transmitter k
% P_HT    : Aggregated RF energy harvesting rate by the device equipped with
%           time switching architecture

% F,F_bar : Used to define thermal noise and noise in channel
%% Initialization and Input parameter settings
K       = 5;
tau     = 0.5;
delta   = 1; %so that 'hk' can be equal to 'exp(theta)'
beta    = 0.3; % conversion efficiency : 30%
mu      = 0.01; %must be grater than 0
c       = zeros(K,1);
dk      = zeros(K,1);
e       = 0.05;
gamma   = 4;
gamma_A = 5;
theta_A = 20;
R       = 30;  % 30 m
Pa      = 39.81; %39.81W =46dBm
Ps      = 10; %100mW
lemda   = 2;
theta   = 0.5; %m/omega
da      = 2;  % 80 m
xie     = [0 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1]; %needs to be variable
%% Input data 
% All data that needed to be entered by the user
%input('Portion of time to harvest energy: ');
%input('RF to DC power conversion efficiency');
%randi([0,MAX],[K,1]); %power value changes from '0' to 'max' with 'K'(no of transmitter)by 'T' (instances of time)
%input('Rate parameter: '); %theta must be positive, for all transmitters(1 by k)
P_H = zeros(K,1);
hk=zeros(K,1);
si=zeros(20,1);
xk =[50 40 70; 60 80 90; 59 92 63; 50 94 133;120 50 31]; %input('Enter coordinates: '); %put array of coordinates for K no of transmitters,'(K,3)'
%% Defining different functions to be used further
%%% Distance of transmitter from device
for l=1:2
    
 for q=1:2
   
     for i=1:K
      c(i) = ((((xk(i,1))^2)+((xk(i,2))^2)+((xk(i,3))^2))^0.5);
      dk(i) = e+c(i);
     end
%%% Channel gain between transmitter and device
    
     for i=1:K
      hk(i) = theta*exp(-theta*dk(i));
     end
%%% Channel gain between Access point and device
     ha  = lemda*exp(lemda*da);
     acc = (Pa*ha)/(da^gamma_A);
%%% Channel and internal Noise approx. function
     F = mu*exp(-mu);
     u = mu*( -mu + 1 ); 
     syms x
     Ei = double(-( int((1/(mu*(x+1)))*exp(-mu*(x+1)),x,mu,inf, 'PrincipalValue',false))); % Calculate value of the integral in the paper mailed page no 9 and conver the result into double 
     F_bar = (1/(-mu*exp(mu)*Ei));

%% Calculating Energy Harvesting rate
%%% RF energy harvesting rate instantaneous
    
     for i=1:K
      P_H(i) = (Ps*hk(i))/((dk(i))^gamma);
     end

%%% Aggregate RF energy harvesting rate
     s = sum(P_H)+acc;
     P_HT= (tau*beta)*((F_bar)/(1+F))*s;  % for the expression before 's' we can have some other expression to show noise

%%% Expectation of RF energy harvesting rate 
     cons = (Pa/((theta_A)*(da^gamma_A)));
     
     for i=1:20
      si(i)=(2*pi*(xie(i))*Ps*delta)/theta;
     end
     
     if (gamma~=1 && gamma~=2 )
       Ir = ((e^(2-gamma))-(((R+e)^(1-gamma))*(e+(gamma-1)*R)))/((gamma-2)*(gamma-1));
     else if (gamma==1)
       Ir = R -(e*log(1+(R/e)));
     else Ir = (log(1+(R/e)))- (R/(R+e));
         end
     end
     
     E_PHTS=zeros(20,1);
     for i=1:20
      E_PHTS(i)= (tau*beta)*(cons+((si(i))*Ir));
     end
%% Results and Plots
%%% Plot for RF energy harvesting rate vs density of Ambient Rf
%%% transmitters
%plot(xie,E_PHTS);
   
     semilogy(xie,E_PHTS,'Linewidth',2);
     hold on
    
   e=e+0.15;
 end
 gamma=gamma+1;
end
