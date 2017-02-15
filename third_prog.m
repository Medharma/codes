%% Program for Optimal Policy for RF energy Harvesting Rate
%%% Network model: Time Switching Architecture
% That is, either the information receiver or the RF energy harvester is connected
% to the antenna at a given time.

%% Variable description
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
% P_H     : RF energy harvesting rate(in watts) by the device from Rf
%           transmitter k
% P_HT    : Aggregated RF energy harvesting rate by the device equipped with
%           time switching architecture

% F,F_bar : Used to define thermal noise and noise in channel
%% Initialization
T     = 10;
K     = 5;
delta = 1; %so that 'hk' can be equal to 'exp(theta)'
max   = 10;
mu    = 3; %must be grater than 0
c     = zeros(K,1);
dk    = zeros(K,1);
e     = [0.1; 0.2; 0.3; 0.4; 0.5]; %different values taken for plots (different scenerios)
gamma = 0.5;  %%[0.1; 0.2; 0.3; 0.4; 0.5];

%% Input data 
% All data that needed to be entered by the user
tau = 0.5;%input('Portion of time to harvest energy: ');
beta = 0.5;%input('RF to DC power conversion efficiency');
Ps = randi([0,max],[K,1]); %power value changes from '0' to 'max' with 'K'(no of transmitter)by 'T' (instances of time)
theta = [21 23 13 54 11];%input('Rate parameter: '); %theta must be positive, for all transmitters(1 by k)

xk = [1 2 3; 1 1 1; 5 2 3; 5 4 3;2 3 1];%input('Enter coordinates: '); %put array of coordinates for K no of transmitters,'(K,3)'
%% Defining different functions to be used further
hk = exp(theta);
F=exp(mu);
syms x;
u=mu*(x+1); %%%%%% expression needs to be verified
f=(1/u)*exp(-u);
Ei=-(int('f',u,-mu,inf)); %, 'IgnoreSpecialCases', true));%not valid for value of mu to be zero: ignore special case can be removed
F_bar=(1/(-mu*exp(mu)*Ei));
for j=1:5  %%%%%$$loop to work for different scenerios
   for i=1:K
   c(j,i) = ((((xk(i,1))^2)+((xk(i,2))^2)+((xk(i,3))^2))^0.5);
   dk(j,i) = e(j,:)+c(j,i); %%%%$$check this one
   end
end
%% Calculating Energy Harvesting rate
%%% RF energy harvesting rate
P_H=[];
for j=1:5
   for i=1:K
   P_H(j,i) = ((tau*beta)*Ps(:,i)*hk(i))/((dk(j,i))^gamma);
   end
end
%%% Aggregate RF energy harvesting rate
s = sum(P_H,2); %%%%%%$$need to sum over same e(sum of each row)
P_HT=((F_bar)/(1+F))*s;  % for the expression before 's' we can have some other expression to show noise

%% Results and plots
% to plot the results of harvesting rate vs density of ambient RF transmitter 
%%% Plotting P_TS vs e
%%% plotting P_TS vs gamma
