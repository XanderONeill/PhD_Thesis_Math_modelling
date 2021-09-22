%%% Deterministic SICS model

function dydt = DETSICS(t, y, alpha)

dydt = zeros(3,1);

N = y(1) + y(2) + y(3);

%% Parameters
b = 10;
d = 1;
K = 1000;
q = (b - d)/(b*K);

beta_F = 1;
beta_D = 0.2;

epsilon = 0;
kappa = 100;

gamma = 50;
%alpha = 50;

%% SIS Model

dydt(1) = b*N*(1-q*N) - beta_F*y(1)*(epsilon*y(3) + y(2))/N - beta_D*y(1)*(epsilon*y(3) + y(2)) - d*y(1) + kappa*y(3);
dydt(2) = beta_F*y(1)*(epsilon*y(3) + y(2))/N + beta_D*y(1)*(epsilon*y(3) + y(2)) - y(2)*(alpha + gamma + d);
dydt(3) = gamma*y(2) - y(3)*(kappa + d);

end
