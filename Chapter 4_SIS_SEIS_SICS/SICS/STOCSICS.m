%%% Stochastic SICS model using Gillespie's Algorithm

function pop = STOCSICS(t, y, alpha, epsilon)
N = sum(y);
S = y(1);
I = y(2);
C = y(3);


%% Parameters
b = 10;
d = 1;
K = 1000;
q = (b - d)/(b*K);
beta_D = 0.15; %0.15;
beta_F = 30; %beta_D*K/5;

kappa = 80;
gamma = 50;

%epsilon = 0.5; 
c=0.5;
alpha_C = c*alpha;
%alpha = 50;

%% MC Run

ev = zeros(7,2);  %event vector
tv = zeros(7,3);  %transition vector

ev1 = max(0, b*N*(1-q*N));
RT = ev1 + beta_F*S*(I + epsilon*C)/N + beta_D*S*(I + epsilon*C) + d*S + kappa*C + I*(alpha + gamma + d) + C*(d+alpha_C);
    
ev(1,:) = [0, ev1/RT];
ev(2,:) = [ev(1,2), ev(1,2) + (beta_F*S*(I + epsilon*C)/N + beta_D*S*(I + epsilon*C))/RT];
ev(3,:) = [ev(2,2), ev(2,2) + d*S/RT];
ev(4,:) = [ev(3,2), ev(3,2) + kappa*C/RT];
ev(5,:) = [ev(4,2), ev(4,2) + gamma*I/RT];
ev(6,:) = [ev(5,2), ev(5,2) + I*(alpha+d)/RT];
ev(7,:) = [ev(6,2), 1];

%set the upper limit to be greater than 1 so that the code is more
%efficient later on in choosing an event. It doens't matter what value this
%takes as all uniform random numbers lie between 0 and 1 inclusive.

tv(1,:) = [ 1, 0, 0];
tv(2,:) = [-1, 1, 0];
tv(3,:) = [-1, 0, 0];
tv(4,:) = [ 1, 0,-1];
tv(5,:) = [ 0,-1, 1];
tv(6,:) = [ 0,-1, 0];
tv(7,:) = [ 0, 0,-1];
 
U = rand(1);
for i = 1:7
    if U >= ev(i,1) && U < ev(i,2)
        y = y + tv(i,:);
        break
    end
end

dt = exprnd(1/RT);
t = t + dt;
 
pop = [t, y];

end