function dydt = Model_tick_host(t, y, H_s, H_u)
dydt = zeros(7,1);
%%% Maximum birth rate
r_T = 2000;                   %Matser et al. 2016

%%% Natural death rates
b_L = 0.0365;                           %1/15;  %Valcarel et al. 2016
b_N = 0.015;                            %1/29;
b_A = 0.00625;                          %1/42;

%%% Rates of Detachment
sigma_L = 0.28;                         %1/10;  %Fran Ruiz Fons
sigma_N = 0.22;                         %1/10;
sigma_A = 0.12;                         %1/20;

%% Rates of attachment

beta_1 = 0.015;                         
beta_2 = 0.0005;                        
beta_3 = 0.03;                          
beta_4 = 0.13;                          

%%% Density dependent terms
%n_S = 95; n_U = 19;
u = 0.04; v = 0.4;
s_T = 0.001/v;                            
s_L = 0.001/u;                           
s_NS = s_L;
s_NU = 0.001/v;                           
    
r_L = 1; r_NS = 1; r_NU = 1;               %maximum probability = 1
a_T = r_T/(1+s_T*(y(7)+y(5))/H_u);
m_L = r_L/(1+s_L*(y(2)+y(4))/H_s);
m_NS= r_NS/(1+s_NS*(y(2)+y(4))/H_s);
m_NU= r_NU/(1+s_NU*(y(7)+y(5))/H_u);

dydt(1) = a_T*sigma_A*y(7) - beta_1*H_s*y(1) - b_L*y(1); 
dydt(2) = beta_1*H_s*y(1)  - sigma_L*y(2);
dydt(3) = m_L*sigma_L*y(2) - (beta_2*H_s + beta_3*H_u)*y(3) - b_N*y(3);
dydt(4) = beta_2*H_s*y(3)  - sigma_N*y(4);
dydt(5) = beta_3*H_u*y(3)  - sigma_N*y(5);
dydt(6) = (m_NS*y(4)+m_NU*y(5))*sigma_N - beta_4*H_u*y(6) - b_A*y(6);
dydt(7) = beta_4*H_u*y(6)  - sigma_A*y(7);
end