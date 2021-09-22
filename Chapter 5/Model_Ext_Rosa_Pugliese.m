function dydt = Model_Ext_Rosa_Pugliese(t, y, H_s, H_u)
dydt = zeros(6,1);

T = y(1) + y(2) + y(3) + y(4) + y(5) + y(6);

r_T = 2000;
b_L = 0.0365;
b_N = 0.015;
b_A = 0.00625;

sigma_A = 0.12;
sigma_L = 0.28;
sigma_N = 0.22;

beta_1 = 0.015;
beta_2 = 0.0005;
beta_3 = 0.03;
beta_4 = 0.13;

r_L = 1;
r_N = 1;
s_T = 0.001;
s_L = s_T;
s_N = s_T;
u = 0.04; v = 0.4;

a_T = r_T/(1 + s_T*T/(u*H_s + v*H_u));
m_L = r_L/(1 + s_L*T/(u*H_s + v*H_u));
m_N = r_N/(1 + s_N*T/(u*H_s + v*H_u));

dydt(1) = sigma_A*a_T*y(6) - beta_1*H_s*y(1) - b_L*y(1);
dydt(2) = beta_1*H_s*y(1) - sigma_L*y(2);
dydt(3) = sigma_L*m_L*y(2) - beta_2*H_s*y(3) - beta_3*H_u*y(3) - b_N*y(3);
dydt(4) = beta_2*H_s*y(3) + beta_3*H_u*y(3) - sigma_N*y(4);
dydt(5) = sigma_N*m_N*y(4) - beta_4*H_u*y(5) - b_A*y(5);
dydt(6) = beta_4*H_u*y(5) - sigma_A*y(6);
end
