function dydt = ModelTBwASFdeath(t, y, prevalence, K1)
dydt = zeros(10,1);

%%% Model function for wild boar TB with added mortality due to ASF.
%%% Function used for analysis section of Chapter 3.

%% The Populations

%First status of infection refers to TB
%Second status of infection refers to ASF

%y(1)->y(3) piglets
%y(4)->y(6) yearlings
%y(7)->y(9) adults

%y(10): free-living TB particles

P = y(1)  + y(2)  + y(3);  %Total Piglets
Y = y(4) + y(5) + y(6); %Total Yearlings
A = y(7) + y(8) + y(9); %Total Adults
N = P + Y + A; %Total Population
G = y(3) + y(6) + y(9); %Total TB Generalised Individuals

%% Parameters

if prevalence == 0
    b = log(4); %birth rate
elseif prevalence == 1
    b = log(7);
end

if prevalence == 0
    K = K1;     
elseif prevalence == 1
    K = K1;   %65.63     
end

m = 1;        %maturity of every other class of pigs

d_P = 1/7;
d_A = d_P;
d_Y = d_A;

q = (1 - (d_A*(d_P+m)*(d_Y+m))/(m*b*(m+d_A)))/K;
%q = (b - d_P)/(b*K);

c = 0.3; %culling rate
%% TB Transmission Dynamics

c_beta = 3;

beta_DA = 1.12;
beta_DP = c_beta*beta_DA;
beta_DY = beta_DP;

beta_FA = 6.5/27;      %25/5.92
beta_FP = c_beta*beta_FA;
beta_FY = beta_FP;

if prevalence == 0 
    omega = 0.1;         %Low Prevalence = 0.1
    c_epsilon = 1;       %Low Prevalence = 1
elseif prevalence == 1
    omega = 1;           %High Prevalence = 1
    c_epsilon = 3;       %High Prevalence = 3
end  
             
epsilon_A = 2/3;                 %rate of transmission from infected to generalised
epsilon_P = c_epsilon*epsilon_A; %rate of transmission from infected to generalised
epsilon_Y = epsilon_P;           %rate of transmission from infected to generalised

alpha = 1; %mortality rate due to TB

%% Free-Living Parameters

lambda = 1; %rate of shedding from generalised class
mu = 6; %rate of decay of free-living particles

%% ASF Death
if prevalence == 0
    ASF = 0.0029;
elseif prevalence == 1
    ASF = 0.0071;
end
gamma = (365*0.85)/5;

%% The Model - Piglets

dydt(1) = b*(Y+A)*(1-q*N) - beta_DP*y(1)*G/N - omega*beta_FP*y(1)*y(10)                  - m*y(1) - d_P*y(1) - ASF*gamma*y(1); 
dydt(2) =                   beta_DP*y(1)*G/N + omega*beta_FP*y(1)*y(10) - epsilon_P*y(2) - m*y(2) - d_P*y(2) - ASF*gamma*y(2);
dydt(3) =                   epsilon_P*y(2) - alpha*y(3) - m*y(3) - d_P*y(3) - ASF*gamma*y(3);

%% Yearlings

dydt(4) = m*y(1) - m*y(4) - d_Y*y(4) - c*y(4) - beta_DY*y(4)*G/N - omega*beta_FY*y(4)*y(10) - ASF*gamma*y(4);                  
dydt(5) = m*y(2) - m*y(5) - d_Y*y(5) - c*y(5) + beta_DY*y(4)*G/N + omega*beta_FY*y(4)*y(10) - epsilon_Y*y(5) - ASF*gamma*y(5);
dydt(6) = m*y(3) - m*y(6) - d_Y*y(6) - c*y(6) + epsilon_Y*y(5) - alpha*y(6) - ASF*gamma*y(6); 

%% Adults

dydt(7) = m*y(4) - d_A*y(7) - c*y(7) - beta_DA*y(7)*G/N - omega*beta_FA*y(7)*y(10) - ASF*gamma*y(7);
dydt(8) = m*y(5) - d_A*y(8) - c*y(8) + beta_DA*y(7)*G/N + omega*beta_FA*y(7)*y(10) - epsilon_A*y(8) - ASF*gamma*y(8);
dydt(9) = m*y(6) - d_A*y(9) - c*y(9) + epsilon_A*y(8) - alpha*y(9) - ASF*gamma*y(9);      

%% The Model - Carcasses and Free-Living

dydt(10) = lambda*G - mu*y(10);

end