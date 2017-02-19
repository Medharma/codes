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
alpha   = -1; %aplha=0 for poisson
sigs    = 10^(-15); %sig^2 = -120dbm/hz = 1*10(-15) W/hz
xie     = 0.005;%[0 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1]; %needs to be variable
h       = zeros(N,1);
d       = zeros(N,1);

%%f_bar
 F = mu*exp(-mu);
     u = mu*( -mu + 1 ); 
     syms v
     Ei = double(-( int((1/(mu*(v+1)))*exp(-mu*(v+1)),v,mu,inf, 'PrincipalValue',true))); % Calculate value of the integral in the paper mailed page no 9 and conver the result into double 
     F_bar = (1/(-mu*exp(mu)*Ei));
%%%x and y
%% Fredholm Determinant using numerical Analysis formula
%for Equation 24
K=xie*(exp(pi*xie*x*y_bar))*(exp((-pi/2)*xie)*(((mod(x))^2)+((mod_y)^2)));
A1=((1-(1+((mu*tau*beta*Ps*F_bar)/(theta*Pc*((d_x+e)^gamma)))))^0.5);
A2=((1-(1+((mu*tau*beta*Ps*F_bar)/(theta*Pc*((d_y+e)^gamma)))))^0.5);
A=A1*K*A2;
Det_24=det((I+alpha*A)^(-1/alpha));

%for Equation 25
Am1=((1-((1+((theta_A*(d_xa^gamma_A)*Ps*Ic*(2^(kappa/(W(1-tau)))-1))/(theta*Pa*Pc*((d_x)+e)^gamma)))^(-delta)))^0.5);  %Ic is interference coefficient
Am2=((1-((1+((theta_A*(d_xa^gamma_A)*Ps*Ic*(2^(kappa/(W(1-tau)))-1))/(theta*Pa*Pc*((d_y)+e)^gamma)))^(-delta)))^0.5);
Am=Am1*K*Am2;

Bm1=((1-((1+(((theta_A*(d_xa^gamma_A)*Ps*Pc*Ic*(2^(kappa/(W(1-tau)))-1))+(mu*tau*beta*Ps*Pa*F_bar))/(theta*Pa*Pc*((d_x)+e)^gamma)))^(-delta)))^0.5);  %Ic is interference coefficient
Bm2=((1-((1+(((theta_A*(d_xa^gamma_A)*Ps*Pc*Ic*(2^(kappa/(W(1-tau)))-1))+(mu*tau*beta*Ps*Pa*F_bar))/(theta*Pa*Pc*((d_y)+e)^gamma)))^(-delta)))^0.5);
Bm=Bm1*K*Bm2;

Det_25_am=det((I+alpha*Am)^(-1/alpha));
Det_25_bm=det((I+alpha*Bm)^(-1/alpha));
