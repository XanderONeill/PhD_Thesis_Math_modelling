function dydt = ASFModelEstoniaINFCHR(t, y, feeding, bF, bE, rho, c)
dydt = zeros(11,1);

%%% function for the ASF model framework representative of Estonia, but
%%% with the survivor class of individuals having the potential to be
%%% infectious

%%% model framework is different to the standard ASFModelEstonia as we look
%%% at the two types of survivors as discussed in Stahl et al. 2019 (see
%%% main paper for citation).

%% The Populations
%y(1) - susceptible piglet
%y(2) - initially infected piglet
%y(3) - "survivor" piglet
%y(4) - secondarily infected piglet
%y(5) - recovered piglet

%y(6) - susceptible adult
%y(7) - initially infected adult
%y(8) - "survivor" adult
%y(9) - secondarily infected adult
%y(10) - recovered adult

%y(11) - infected carcasses

% S -> II -> C -> I
%              -> R

P = y(1) + y(2) + y(3) + y(4) + y(5);
A = y(6) + y(7) + y(8) + y(9) + y(10);
N = P + A;                             %total alive population
I = y(2) + y(4) + y(7) + y(9);
C = y(3) + y(8);


%% Basic Parameters

%Carrying Capacities
if feeding == 0
    K = 2;     %Estonia = 2, Spain = 5
elseif feeding == 1
    K = 4;     %Estonia = 4, Spain = 10
end

%average value of birth rate a(t)
if feeding == 0
    avga = 3;
elseif feeding == 1
    avga = 6;
end

% The birth rate a(t)
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

alpha = 12/10;  %rate of maturity

b_A = 0.33;         %adult natural death rate
b_H = 0.49;         %adult hunting mortality rate
b_D = b_A + b_H;
b_0 = log(2);         %piglet natural death rate at low densities

%additional density dependent term for piglet natural death rate
if feeding == 0
    b_1 = 1/1.33*((avga - b_D)*alpha - b_0*(b_D))/(b_D*K);
elseif feeding == 1
    b_1 = 1/1.08695*((avga - b_D)*alpha - b_0*(b_D))/(b_D*K);
end

b_P = (b_0 + b_1*N);     %total piglet natural death rate

%% Control Measures

b_C = 0;       %culling rate
r = 0;         %carcass removal rate

%% Transmission Parameters

beta_F = bF;     %Direct contact transmission - freq. dep.
%c = 0.3;         %Proportion of transmission from survivor individuals
beta_E = bE;      %Environmental transmission - dens. dep.
gamma = 365/5;   %disease induced mortality rate
                 %5 days
%rho = 0.15;      %proportion progressing to the "survivor" class
eta = 9.125;     %rate of progression to recovered class
                 %30-40 days
kappa = 12/6;    %rate at which "survivors" progress to further infected class
                 %6 months
%Can replace kappa with (1-phi)*eta/phi where phi represents the proportion
%of "survivors" which will progress to recovered class 0 < phi < 1

mu_C = 52/8;     %Carcass degradation rate
                 %Estonia - 8 weeks approx

%% The Model

%Piglets
dydt(1) = a*A - beta_F*y(1)*(I+c*C)/N - beta_E*y(1)*y(11) - y(1)*(alpha + b_P + b_C);
dydt(2) =       beta_F*y(1)*(I+c*C)/N + beta_E*y(1)*y(11) - y(2)*(alpha + b_P + b_C + gamma);
dydt(3) =       gamma*rho*y(2)                            - y(3)*(alpha + b_P + b_C + eta + kappa);   
dydt(4) =       kappa*y(3)                                - y(4)*(alpha + b_P + b_C + gamma);
dydt(5) =       eta*y(3)                                  - y(5)*(alpha + b_P + b_C);

%Yearlings and Adults
dydt(6) =  alpha*y(1) - beta_F*y(6)*(I+c*C)/N - beta_E*y(6)*y(11) - y(6)*(b_A + b_H + b_C);
dydt(7) =  alpha*y(2) + beta_F*y(6)*(I+c*C)/N + beta_E*y(6)*y(11) - y(7)*(b_A + b_H + b_C + gamma);
dydt(8) =  alpha*y(3) + gamma*rho*y(7)                            - y(8)*(b_A + b_H + b_C + eta + kappa);   %(1-phi)*eta/phi
dydt(9) =  alpha*y(4) + kappa*y(8)                                - y(9)*(b_A + b_H + b_C + gamma);
dydt(10) = alpha*y(5) + eta*y(8)                                  - y(10)*(b_A + b_H + b_C);

%The Dead
dydt(11) = gamma*(1-rho)*(y(2)+y(7)) + gamma*(y(4)+y(9)) + b_P*(y(2)+y(3)+y(4)+y(5)) + b_A*(y(7)+y(8)+y(9)+y(10)) - mu_C*y(11) - r*y(11);
end