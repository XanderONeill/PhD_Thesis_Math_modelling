function dydt = ModelTBandASF(t, y, prevalence, d1, b_C, r)
dydt = zeros(29,1);

%%% Model function representative of the TB and ASF system in Spain.
%%% Takes in the following inputs:

%prevalence = 0 or 1 
    %where 0 denotes a low prevalence region and 1 denotes a high prevalence region
%d1 - degradation rate of carcasses (per year)
%b_C is the strength of culling
    %if set to 0, then no culling is present
%r is the carcass removal strength
    %if set to 0, then no carcass removal is present

%% The Populations

%y(1)->y(9) piglets
%y(10)->y(18) yearlings
%y(19)->y(27) adults

%y(1): susceptible, susceptible   y(2): infected, susceptible  
%y(3): generalised, susceptible   y(4): susceptible, infected      
%y(5): infected, infected         y(6): generalised, infected   
%y(7): susceptible, chronic       y(8): infected, chronic  
%y(9): generalised, chronic 

%y(10): susceptible, susceptible   y(11): infected, susceptible  
%y(12): generalised, susceptible   y(13): susceptible, infected      
%y(14): infected, infected         y(15): generalised, infected   
%y(16): susceptible, chronic       y(17): infected, chronic  
%y(18): generalised, chronic   

%y(19): susceptible, susceptible   y(20): infected, susceptible  
%y(21): generalised, susceptible   y(22): susceptible, infected      
%y(23): infected, infected         y(24): generalised, infected   
%y(25): susceptible, chronic       y(26): infected, chronic  
%y(27): generalised, chronic 

%y(28): infected carcasses
%y(29): free-living TB particles

P = y(1)  + y(2)  + y(3)  + y(4)  + y(5)  + y(6)  + y(7)  + y(8)  + y(9);  %Total Piglets
Y = y(10) + y(11) + y(12) + y(13) + y(14) + y(15) + y(16) + y(17) + y(18); %Total Yearlings
A = y(19) + y(20) + y(21) + y(22) + y(23) + y(24) + y(25) + y(26) + y(27); %Total Adults
N = P + Y + A; %Total Population

I = y(4) + y(5) + y(6) + y(13) + y(14) + y(15) + y(22) + y(23) + y(24); %Total ASF Infected Individuals
G = y(3) + y(6) + y(9) + y(12) + y(15) + y(18) + y(21) + y(24) + y(27); %Total TB Generalised Individuals

%% Parameters

if prevalence == 0
    b = log(4); %birth rate
elseif prevalence == 1
    b = log(7);
end

if prevalence == 0
    K = 5.95;     
elseif prevalence == 1
    K = 27;   %65.63     
end

m = 1;        %maturity of every other class of pigs

d_P = 1/7;
d_A = d_P;
d_Y = d_A;

q = (1 - (d_A*(d_P+m)*(d_Y+m))/(m*(b*m+b*d_A)))/K;

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

%% ASF Transmission Dynamics

beta_F = 56; %frequency dependent transmission rate
beta_E = 2.9; %environmental transmission rate
gamma = 365/5; %mortality rate due to ASF
kappa = 2; %rate of re-infection from chronic class
rho = 0.85; %proportion of those dying from ASF
%(1-rho) - proportion of infected that don't die but enter the chronic
%class

%% Carcass Paramters

d = d1; %degradation rate of carcasses
%Low Degradation Rate = 52/4
%High Degradation Rate = 52/1

%% Free-Living Parameters

lambda = 1; %rate of shedding from generalised class
mu = 6; %rate of decay of free-living particles

%% The Model
    
%Piglets 
dydt(1) = b*(Y+A)*(1-q*N) - beta_DP*y(1)*G/N - omega*beta_FP*y(1)*y(29)                  - m*y(1) - d_P*y(1) - beta_F*y(1)*I/N - beta_E*y(1)*y(28) - b_C*y(1);
dydt(2) =                   beta_DP*y(1)*G/N + omega*beta_FP*y(1)*y(29) - epsilon_P*y(2) - m*y(2) - d_P*y(2) - beta_F*y(2)*I/N - beta_E*y(2)*y(28) - b_C*y(2);
dydt(3) =                   epsilon_P*y(2) - alpha*y(3)                                  - m*y(3) - d_P*y(3) - beta_F*y(3)*I/N - beta_E*y(3)*y(28) - b_C*y(3);

dydt(4) = beta_F*y(1)*I/N + beta_E*y(1)*y(28) - m*y(4) - d_P*y(4) - beta_DP*y(4)*G/N - omega*beta_FP*y(4)*y(29)                  - gamma*y(4) + kappa*y(7) - b_C*y(4);
dydt(5) = beta_F*y(2)*I/N + beta_E*y(2)*y(28) - m*y(5) - d_P*y(5) + beta_DP*y(4)*G/N + omega*beta_FP*y(4)*y(29) - epsilon_P*y(5) - gamma*y(5) + kappa*y(8) - b_C*y(5);
dydt(6) = beta_F*y(3)*I/N + beta_E*y(3)*y(28) - m*y(6) - d_P*y(6) + epsilon_P*y(5) - alpha*y(6)                                  - gamma*y(6) + kappa*y(9) - b_C*y(6);

dydt(7) = gamma*(1-rho)*y(4) - kappa*y(7) - m*y(7) - d_P*y(7) - beta_DP*y(7)*G/N - omega*beta_FP*y(7)*y(29)                  - b_C*y(7);
dydt(8) = gamma*(1-rho)*y(5) - kappa*y(8) - m*y(8) - d_P*y(8) + beta_DP*y(7)*G/N + omega*beta_FP*y(7)*y(29) - epsilon_P*y(8) - b_C*y(8);
dydt(9) = gamma*(1-rho)*y(6) - kappa*y(9) - m*y(9) - d_P*y(9) + epsilon_P*y(8) - alpha*y(9)                                  - b_C*y(9);

%Yearlings
dydt(10) = m*y(1) - m*y(10) - d_Y*y(10) - c*y(10) - beta_DY*y(10)*G/N - omega*beta_FY*y(10)*y(29)                   - beta_F*y(10)*I/N - beta_E*y(10)*y(28) - b_C*y(10);
dydt(11) = m*y(2) - m*y(11) - d_Y*y(11) - c*y(11) + beta_DY*y(10)*G/N + omega*beta_FY*y(10)*y(29) - epsilon_Y*y(11) - beta_F*y(11)*I/N - beta_E*y(11)*y(28) - b_C*y(11);
dydt(12) = m*y(3) - m*y(12) - d_Y*y(12) - c*y(12) + epsilon_Y*y(11) - alpha*y(12)                                   - beta_F*y(12)*I/N - beta_E*y(12)*y(28) - b_C*y(12);

dydt(13) = m*y(4) - m*y(13) - d_Y*y(13) - c*y(13) - beta_DY*y(13)*G/N - omega*beta_FY*y(13)*y(29)                   + beta_F*y(10)*I/N + beta_E*y(10)*y(28) - gamma*y(13) + kappa*y(16) - b_C*y(13);
dydt(14) = m*y(5) - m*y(14) - d_Y*y(14) - c*y(14) + beta_DY*y(13)*G/N + omega*beta_FY*y(13)*y(29) - epsilon_Y*y(14) + beta_F*y(11)*I/N + beta_E*y(11)*y(28) - gamma*y(14) + kappa*y(17) - b_C*y(14);
dydt(15) = m*y(6) - m*y(15) - d_Y*y(15) - c*y(15) + epsilon_Y*y(14) - alpha*y(15)                                   + beta_F*y(12)*I/N + beta_E*y(12)*y(28) - gamma*y(15) + kappa*y(18) - b_C*y(15);

dydt(16) = m*y(7) - m*y(16) - d_Y*y(16) - c*y(16) - beta_DY*y(16)*G/N - omega*beta_FY*y(16)*y(29)                   + gamma*(1-rho)*y(13) - kappa*y(16) - b_C*y(16);
dydt(17) = m*y(8) - m*y(17) - d_Y*y(17) - c*y(17) + beta_DY*y(16)*G/N + omega*beta_FY*y(16)*y(29) - epsilon_Y*y(17) + gamma*(1-rho)*y(14) - kappa*y(17) - b_C*y(17);
dydt(18) = m*y(9) - m*y(18) - d_Y*y(18) - c*y(18) + epsilon_Y*y(17) - alpha*y(18)                                   + gamma*(1-rho)*y(15) - kappa*y(18) - b_C*y(18);

%Adults
dydt(19) = m*y(10) - d_A*y(19) - c*y(19) - beta_DA*y(19)*G/N - omega*beta_FA*y(19)*y(29)                   - beta_F*y(19)*I/N - beta_E*y(19)*y(28) - b_C*y(19);
dydt(20) = m*y(11) - d_A*y(20) - c*y(20) + beta_DA*y(19)*G/N + omega*beta_FA*y(19)*y(29) - epsilon_A*y(20) - beta_F*y(20)*I/N - beta_E*y(20)*y(28) - b_C*y(20);
dydt(21) = m*y(12) - d_A*y(21) - c*y(21) + epsilon_A*y(20) - alpha*y(21)                                   - beta_F*y(21)*I/N - beta_E*y(21)*y(28) - b_C*y(21);

dydt(22) = m*y(13) - d_A*y(22) - c*y(22) - beta_DA*y(22)*G/N - omega*beta_FA*y(22)*y(29)                   + beta_F*y(19)*I/N + beta_E*y(19)*y(28) - gamma*y(22) + kappa*y(25) - b_C*y(22);
dydt(23) = m*y(14) - d_A*y(23) - c*y(23) + beta_DA*y(22)*G/N + omega*beta_FA*y(22)*y(29) - epsilon_A*y(23) + beta_F*y(20)*I/N + beta_E*y(20)*y(28) - gamma*y(23) + kappa*y(26) - b_C*y(23);
dydt(24) = m*y(15) - d_A*y(24) - c*y(24) + epsilon_A*y(23) - alpha*y(24)                                   + beta_F*y(21)*I/N + beta_E*y(21)*y(28) - gamma*y(24) + kappa*y(27) - b_C*y(24);

dydt(25) = m*y(16) - d_A*y(25) - c*y(25) - beta_DA*y(25)*G/N - omega*beta_FA*y(25)*y(29)                   + gamma*(1-rho)*y(22) - kappa*y(25) - b_C*y(25);
dydt(26) = m*y(17) - d_A*y(26) - c*y(26) + beta_DA*y(25)*G/N + omega*beta_FA*y(25)*y(29) - epsilon_A*y(26) + gamma*(1-rho)*y(23) - kappa*y(26) - b_C*y(26);
dydt(27) = m*y(18) - d_A*y(27) - c*y(27) + epsilon_A*y(26) - alpha*y(27)                                   + gamma*(1-rho)*y(24) - kappa*y(27) - b_C*y(27);

I_P = y(4) + y(5) + y(6);
I_Y = y(13) + y(14) + y(15);
I_A = y(22) + y(23) + y(24);

%Free-Living and Carcasses
dydt(28) = d_P*I_P + d_Y*I_Y + d_A*I_A + gamma*rho*(I_P + I_Y + I_A) - d*y(28) - r*y(28);
dydt(29) = lambda*G - mu*y(29); 

end