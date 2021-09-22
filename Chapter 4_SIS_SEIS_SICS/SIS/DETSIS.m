function dydt = DETSIS(t, y, alpha, gamma)

dydt = zeros(2,1);

N = y(1) + y(2);

%% Parameters
b = 10;
d = 1;
K = 1000;
q = (b - d)/(b*K);
beta_D = 0.15;
beta_F = beta_D*K/5;

%gamma = 10;
%alpha = 50;

%% SIS Model

dydt(1) = b*N*(1-q*N) - beta_F*y(1)*y(2)/N - beta_D*y(1)*y(2) - d*y(1) + gamma*y(2);
dydt(2) = beta_F*y(1)*y(2)/N + beta_D*y(1)*y(2) - y(2)*(gamma + alpha + d);

end
