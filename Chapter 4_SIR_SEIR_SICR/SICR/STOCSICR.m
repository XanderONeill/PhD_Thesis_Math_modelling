%%% Stochastic SICR model using Gillespie's Algorithm

function pop = STOCSICR(t, y, alpha, epsilon)
N = sum(y);
S = y(1);
I = y(2);
C = y(3);
R = y(4);


%% Parameters
b = 10;
d = 1;
K = 1000;
q = (b - d)/(b*K);
beta_D = 0.15;
beta_F = 30;

kappa = 80;
gamma = 50;
%alpha = 50;

%epsilon = 0;
c = 0.5; 
alpha_C = c*alpha;

%% MC Run

ev = zeros(8,2);  %event vector
tv = zeros(8,4);  %transition vector

ev1 = max(0, b*N*(1-q*N));
RT = ev1 + beta_F*S*(I+epsilon*C)/N + beta_D*S*(I+epsilon*C) + d*S + I*(alpha + gamma + d) + C*(alpha_C + kappa + d) + R*d;
    
ev(1,:) = [0, ev1/RT];
ev(2,:) = [ev(1,2), ev(1,2) + d*S/RT];
ev(3,:) = [ev(2,2), ev(2,2) + (beta_F*S*(I+epsilon*C)/N + beta_D*S*(I+epsilon*C))/RT];
ev(4,:) = [ev(3,2), ev(3,2) + I*(alpha+d)/RT];
ev(5,:) = [ev(4,2), ev(4,2) + gamma*I/RT];
ev(6,:) = [ev(5,2), ev(5,2) + (alpha_C + d)*C/RT];
ev(7,:) = [ev(6,2), ev(6,2) + kappa*C/RT];
ev(8,:) = [ev(7,2), 1];

%set the upper limit to be greater than 1 so that the code is more
%efficient later on in choosing an event. It doens't matter what value this
%takes as all uniform random numbers lie between 0 and 1 inclusive.

tv(1,:) = [ 1, 0, 0, 0];
tv(2,:) = [-1, 0, 0, 0];
tv(3,:) = [-1, 1, 0, 0];
tv(4,:) = [ 0,-1, 0, 0];
tv(5,:) = [ 0,-1, 1, 0];
tv(6,:) = [ 0, 0,-1, 0];
tv(7,:) = [ 0, 0,-1, 1];
tv(8,:) = [ 0, 0, 0,-1];
 
U = rand(1);
for i = 1:8
    if U >= ev(i,1) && U < ev(i,2)
        y = y + tv(i,:);
        break
    end
end

dt = exprnd(1/RT);
t = t + dt;
 
pop = [t, y];

end