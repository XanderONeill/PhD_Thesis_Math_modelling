function dydt = ASFModelSpain(t, y, feeding)
dydt = zeros(8,1);

%feeding = 0 represents natural conditions
%feeding = 1 represents supplemented feeding

%y(1): susceptible piglets  
%y(2): susceptible adults
%y(3): infected piglets  
%y(4): infected adults
%y(5): chronic piglets
%y(6): chronic adults
%y(7): dead (infected) carcasses
%y(8): dead (non-infected) carcasses

A = y(2) + y(4) + y(6); %Adult Population
I = y(3) + y(4); %Infected Population
N = y(1) + y(2) + y(3) + y(4) + y(5) + y(6); %Total alive population

%% Wild Boar Parameters:

%Alpha: Rate of change to adulthood
%Gamma: Mortality/recovery rate;
%Kappa: Rate of reinfection;
%Rho: proportion of infected that die due to the disease
%d: rate of degredation of the carcasses
%K: Estonia Carrying Capacity

alpha = 12/10;
gamma = 365/5; 
kappa= 12/6;       %base line 12/6 
rho = 0.85;        %base line 0.85

%Fixed Deaths
b_D = log(1/0.45);
b_A = 4/10*b_D;
b_0 = log(2);

%Carrying Capacities
% Natural conditions: K=5,  Supplemented feeding: K=10
if feeding == 0
    K = 5;     
elseif feeding == 1
    K = 10;     
end

%average value of birth rate a(t)
if feeding == 0
    avga = 3;
elseif feeding == 1
    avga = 6;
end

%% Control Parameters
b_C = 0;
r = 0;

%% birth

% the seasonal birth function given in Appendix A of the thesis
T = t - floor(t);
if feeding == 0
    scale0 = 3/0.253085;
    a = scale0*exp(-49*(T-(4/12))^2);
elseif feeding == 1
    scale1 = 6/0.4723295;
    if 0 <= T && T < 3/12
        a = scale1*exp(-100*(T-3/12)^2);
    elseif 3/12 <= T && T < 6/12
        a = scale1*(0.5*exp(-106.0903549*(T-3/12)^2) + 0.5);
    elseif 6/12 <= T && T < 9/12
        a = scale1*(0.25*exp(-95*(T-9/12)^2) + 0.5);
    else
        a = scale1*0.75*exp(-95.39708684*(T-9/12)^2);
    end
end

% feeding adjusted for natural conditions/supplemented feeding
if feeding == 0
    b_1 = 1/1.33*((avga - b_D)*alpha - b_0*(b_D))/(b_D*K);
elseif feeding == 1
    b_1 = 1/1.08695*((avga - b_D)*alpha - b_0*(b_D))/(b_D*K);
end

b_P = (b_0 + b_1*N);
 
%% Transmission terms

%beta_F - frequency dependent transmission
%beta_E - density dependent environmental

beta_F = 63;
beta_E = 2;

%In Spain no over-wintered carcasses but quicker degradation
d = 52/1;       %Base Line once per week in Spain

%% The Model

dydt(1) = a*A - beta_F*(y(1)/N)*I - beta_E*y(1)*y(7) - alpha*y(1) - b_P*y(1) - b_C*y(1);
dydt(2) =     - beta_F*(y(2)/N)*I - beta_E*y(2)*y(7) + alpha*y(1) - b_D*y(2) - b_C*y(2);

dydt(3) = beta_F*(y(1)/N)*I + beta_E*y(1)*y(7) - alpha*y(3) - gamma*y(3) + kappa*y(5) - b_P*y(3) - b_C*y(3);
dydt(4) = beta_F*(y(2)/N)*I + beta_E*y(2)*y(7) + alpha*y(3) - gamma*y(4) + kappa*y(6) - b_D*y(4) - b_C*y(4);

dydt(5) = gamma*(1-rho)*y(3) - kappa*y(5) - alpha*y(5) - b_P*y(5) - b_C*y(5);
dydt(6) = gamma*(1-rho)*y(4) - kappa*y(6) + alpha*y(5) - b_D*y(6) - b_C*y(6);

dydt(7) = gamma*rho*I + b_P*y(3) + b_A*y(4) - d*y(7) - r*y(7);
dydt(8) = b_P*y(1) + b_A*(y(2)+y(6)) + b_P*y(5) - r*y(8);

end