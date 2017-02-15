%%Program for MDP;
%% T(theta,a) : immediate throughput function
%  na         : Probability of selected channel to be idle
%  sigma_a    : Probability of successful transmission by secondary user on selected channel a
%  e          : Energy level of energy queue
%  W          : Energy limit
%  q          : Number of packets in data queue

%% starting with assigning probabilities
t=10;
% initilizing 
E_T    = zeros(t,1);
T      = zeros(t,10);
na     = zeros(t,10);
sigma_a = zeros(t,10);

for j = 1:10
    na(j,:) = sort( rand(1,10));
end
for j = 1:10
    sigma_a(j,:)= rand(1, 10);
end 

%% defining Queue state
W = 50;
e = input('Energy level of energy state:  ');
q = input('Number of packets in data queue:  ');
%% Computing immediate throughput
% for different instances
for i = 1:t
    if((e>=W)&&(q>0))
        T(i,:) = na(i,:).*sigma_a(i,:);
    else
        T(i,:) = zeros(1,10);
    end
end
%% Expectation of immediate throughputs

for i=1:t
E_T(i,1) = mean(T(i,:));
end
%% Final expression
fun = sum(E_T)/t;
