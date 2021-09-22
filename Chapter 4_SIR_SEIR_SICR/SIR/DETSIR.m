%%% Deterministic SIR model

function dydt = DETSIR(t, y, alpha)

dydt = zeros(3,1);

N = y(1) + y(2) + y(3);

%% Parameters
b = 10;
d = 1;
K = 1000;
q = (b - d)/(b*K);

beta_F = 1;
beta_D = 0.2;

gamma = 10;
%alpha = 1;

%% SIS Model

dydt(1) = b*N*(1-q*N) - beta_F*y(1)*y(2)/N - beta_D*y(1)*y(2) - d*y(1);
dydt(2) = beta_F*y(1)*y(2)/N + beta_D*y(1)*y(2) - y(2)*(alpha + gamma + d);
dydt(3) = gamma*y(2) - y(3)*d;

end
