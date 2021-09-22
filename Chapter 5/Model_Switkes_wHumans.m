function dydt = Model_Switkes_wHumans(t, y, K, n) %K = 2..100 default K=10km-2
dydt = zeros(8,1);

H = y(1) + y(2) + y(3);
V = y(4) + y(5);
Hu = y(6) + y(7) + y(8);

%K = 10; 
%n = 50; %good quality
%n = 100; %poor quality
aV = 150;   bV = 1;   qV = (aV - bV)/(n*H);
aH = log(4); bH = 1/7; qH = (aH - bH)/K;

betaHV = 0.0112;
betaVV = betaHV/10; betaVH = 50*betaHV;
gammaH = 0.4;
p = 0.1;

aHu = 1/50; betaVHu = betaVV;

dydt(1) = (aH - qH*H)*H - bH*y(1) - betaHV*y(1)*y(5)/H;
dydt(2) = betaHV*y(1)*y(5)/H - bH*y(2) - gammaH*y(2);
dydt(3) = gammaH*y(2) - bH*y(3);

dydt(4) = (aV - qV*V)*(V - p*y(5)) - bV*y(4) - betaVV*y(4)*y(5)/V - betaVH*y(2)*y(4)/H;
dydt(5) = (aV - qV*V)*p*y(5) - bV*y(5) + betaVV*y(4)*y(5)/V + betaVH*y(2)*y(4)/H;

dydt(6) = aHu*Hu - betaVHu*y(6)*y(5)/Hu - aHu*y(6);
dydt(7) = betaVHu*y(6)*y(5)/Hu - gammaH*y(7) - aHu*y(7);
dydt(8) = gammaH*y(7) - aHu*y(8);
end