function dydt = DETSEIS(t, y, alpha, beta_F)
%%% deterministic SEIS model, parameters given in thesis, Chapter 4

dydt = zeros(3,1);

N = y(1) + y(2) + y(3);

%% Natual Dynamics
b = 10;
d = 1;
K = 1000;
q = (b - d)/(b*K);


%% Infection Dynamics
%beta_F = 1;
beta_D = 0.2;

epsilon = 0;
kappa = 100;

gamma = 50;
%alpha = 50;


%% SIS Model

dydt(1) = b*N*(1-q*N) - beta_F*y(1)*(epsilon*y(2) + y(3))/N - beta_D*y(1)*(epsilon*y(2) + y(3)) - d*y(1) + gamma*y(3);
dydt(2) = beta_F*y(1)*(epsilon*y(2) + y(3))/N + beta_D*y(1)*(epsilon*y(2) + y(3)) - y(2)*(kappa + d);
dydt(3) = kappa*y(2) - y(3)*(alpha + gamma + d);

end
