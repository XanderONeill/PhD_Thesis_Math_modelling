function dydt = Model_Modified_Switkes(t, y, K, n) %K = 2..100 default K=10km-2

%%% Modified Switkes model in the absence of humans
%%% Parameters as discussed in Chapter 5 of the thesis

dydt = zeros(5,1);

H = y(1) + y(2) + y(3);
V = y(4) + y(5);

%K = 10;
%n = 50; %good quality
%n = 100; %poor quality
aV = 150;   bV = 1;   qV = (aV - bV)/(n*H);
aH = log(4); bH = 1/7; qH = (aH - bH)/(K);

betaHV = 0.0112;
betaVV = betaHV/10; betaVH = 50*betaHV;
gammaH = 0.4;
p = 0.1;

% Rt = p + betaVV/bV;
% Rl = (betaHV*betaVH*(V/H))/((bH + gammaH)*bV);
% Rtl = Rt/2 + sqrt((Rt^2)/4 + Rl);

dydt(1) = (aH - qH*H)*H - betaHV*y(1)*y(5)/H - bH*y(1);
dydt(2) = betaHV*y(1)*y(5)/H   - gammaH*y(2) - bH*y(2);
dydt(3) =                        gammaH*y(2) - bH*y(3);

dydt(4) = (aV - qV*V)*(V - p*y(5)) - betaVV*y(4)*y(5)/V - betaVH*y(2)*y(4)/H - bV*y(4);
dydt(5) = (aV - qV*V)*p*y(5)       + betaVV*y(4)*y(5)/V + betaVH*y(2)*y(4)/H - bV*y(5);
end