function dydt = Model_tick_host_wInf(t, y, K_S, K_U)
dydt = zeros(20,1);

%%% Populations
H_s = y(15)+y(16)+y(17); %total small mammal density
H_u = y(18)+y(19)+y(20); %total ungulate density

I_TS = y(4) + y(8);                  %Infected ticks on small mammals
T_TS = y(3) + y(4) + y(7) + y(8);    %total ticks on small mammals

I_TU = y(10) + y(14);                %Infected ticks on ungulates
T_TU = y(9) + y(10) + y(13) + y(14); %total ticks on ungulates

%% Tick Non-Infection Parameters
%%% Maximum birth rate
r_T = 2000;                   %Matser et al. 2016

%%% Natural death rates
b_L = 0.0365;                           %1/15;  %Valcarel et al. 2016
b_N = 0.015;                            %1/29;
b_A = 0.00625;  

%%% Rates of Detachment
sigma_L = 0.28;                         %1/10;  %Fran Ruiz Fons
sigma_N = 0.22;                         %1/10;
sigma_A = 0.12;                         %1/20;

%%% Rates of attachment

beta_1 = 0.015;                         %0.00005;
beta_2 = 0.0005;                        %0.0005;
beta_3 = 0.03;                          %0.00008;
beta_4 = 0.13;                          %0.000142;

%%% Density dependent terms
%n_S = 95; n_U = 19;
u = 0.04; v = 0.4;
s_T = 0.001/v;                            %1.996474736/n_U; %Rosa and Pugliese 2007
s_L = 0.001/u;                            %0.5431756071/n_S; 
s_NS = s_L;
s_NU = 0.001/v;                           %0.6598754496/n_U;
    
r_L = 1; r_NS = 1; r_NU = 1;               %maximum probability = 1
a_T = r_T/(1+s_T*T_TU/H_u);
m_L = r_L/(1+s_L*T_TS/H_s);
m_NS= r_NS/(1+s_NS*T_TS/H_s);
m_NU= r_NU/(1+s_NU*T_TU/H_u);

%% Hosts Non-Infection Parameters
a_S = 0.1;
a_U = 0.05;
b_S = 1/1825;
b_U = 1/5475;
q_S = (a_S - b_S)/(a_S*K_S);
q_U = (a_U - b_U)/(a_U*K_U);

%% Infection Parameters
rho = 0.2; %vertical transmission

pL = 0.2;   %infected small mammal -> larvae transmission probability
pN1 = 0.2;  %infected small mammal -> nymph transmission probability
pN2 = 0.2;  %infected ungulate -> nymph transmission probability
pA = 0.2;   %infected ungulate -> adult transmission probability
qL =  0.000005;  %infected larvae -> small mammal transmission probability
qN1 = 0.000005; %infected nymph -> small mammal transmission probability
qN2 = 0.0001; %infected nymph -> ungulate transmission probability
qA =  0.0001;  %infected adult -> ungulate transmission probability

theta_TT = 0.1; %tick -> tick cofeeding transmission coefficient

gamma_S = 0.01; %small mammal recovery rate
gamma_U = 0.01; %ungulate recovery rate


%% The Model

dydt(1)  = a_T*sigma_A*(y(13) + (1-rho)*y(14)) - beta_1*H_s*y(1) - b_L*y(1); 
dydt(2)  = a_T*sigma_A*(rho*y(14))             - beta_1*H_s*y(2) - b_L*y(2);
dydt(3)  = beta_1*(H_s - pL*y(16))*y(1)           - theta_TT*I_TS/T_TS*y(3)  - sigma_L*y(3);
dydt(4)  = beta_1*H_s*y(2) + pL*beta_1*y(16)*y(1) + theta_TT*I_TS/T_TS*y(3)  - sigma_L*y(4);
dydt(5)  = m_L*sigma_L*y(3)  - (beta_2*H_s + beta_3*H_u)*y(5) - b_N*y(5);
dydt(6)  = m_L*sigma_L*y(4)  - (beta_2*H_s + beta_3*H_u)*y(6) - b_N*y(6);
dydt(7)  = beta_2*(H_s - pN1*y(16))*y(5)           - theta_TT*y(7)*I_TS/T_TS  - sigma_N*y(7);
dydt(8)  = beta_2*H_s*y(6) + pN1*beta_2*y(16)*y(5) + theta_TT*y(7)*I_TS/T_TS  - sigma_N*y(8);
dydt(9)  = beta_3*(H_u - pN2*y(19))*y(5)           - theta_TT*y(9)*I_TU/T_TU  - sigma_N*y(9);
dydt(10) = beta_3*H_u*y(6) + pN2*beta_3*y(19)*y(5) + theta_TT*y(9)*I_TU/T_TU  - sigma_N*y(10);
dydt(11) = sigma_N*(m_NS*y(7) + m_NU*y(9))  - beta_4*H_u*y(11) - b_A*y(11);
dydt(12) = sigma_N*(m_NS*y(8) + m_NU*y(10)) - beta_4*H_u*y(12) - b_A*y(12);
dydt(13) = beta_4*(H_u - pA*y(19))*y(11)            - theta_TT*I_TU/T_TU  - sigma_A*y(13);
dydt(14) = beta_4*H_u*y(12) + pA*beta_4*y(19)*y(11) + theta_TT*I_TU/T_TU  - sigma_A*y(14);

dydt(15) = a_S*H_s*(1 - q_S*H_s) - (qL*beta_1*y(2) + qN1*beta_2*y(6))*y(15)     - b_S*y(15);
dydt(16) = (qL*beta_1*y(2) + qN1*beta_2*y(6))*y(15)     - gamma_S*y(16)        - b_S*y(16);
dydt(17) =                                              gamma_S*y(16)        - b_S*y(17);

dydt(18) = a_U*H_u*(1 - q_U*H_u) - (qN2*beta_3*y(6) + qA*beta_4*y(12))*y(18)    - b_U*y(18);
dydt(19) = (qN2*beta_3*y(6) + qA*beta_4*y(12))*y(18)     - gamma_U*y(19)        - b_U*y(19);
dydt(20) =                                              gamma_U*y(19)        - b_U*y(20);
end