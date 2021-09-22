%%% Deterministic SEIR model

function dydt = DETSEIR(t, y, alpha)

dydt = zeros(4,1);

N = y(1) + y(2) + y(3) + y(4);

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
%alpha = 1;

%% SIS Model

dydt(1) = b*N*(1-q*N) - beta_F*y(1)*(epsilon*y(2) + y(3))/N - beta_D*y(1)*(epsilon*y(2) + y(3)) - d*y(1);
dydt(2) = beta_F*y(1)*(epsilon*y(2) + y(3))/N + beta_D*y(1)*(epsilon*y(2) + y(3)) - y(2)*(kappa + d);
dydt(3) = kappa*y(2) - y(3)*(alpha + gamma + d);
dydt(4) = gamma*y(3) - d*y(4);

end
