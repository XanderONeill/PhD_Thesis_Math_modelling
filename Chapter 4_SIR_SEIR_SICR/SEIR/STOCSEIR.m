%%% Stochastic SEIR model using Gillespie's Algorithm

function pop = STOCSEIR(t, y, alpha, epsilon)
N = sum(y);
S = y(1);
E = y(2);
I = y(3);
R = y(4);


%% Parameters
b = 10;
d = 1;
K = 1000;
q = (b - d)/(b*K);
beta_D = 0.15;
beta_F = 30;

%epsilon = 0;
kappa = 80;

gamma = 50;
%alpha = 50;

%% MC Run

ev = zeros(8,2);  %event vector
tv = zeros(8,4);  %transition vector

ev1 = max(0, b*N*(1-q*N));
RT = ev1 + beta_F*S*(I+epsilon*E)/N + beta_D*S*(I+epsilon*E) + d*S + E*(kappa + d) + I*(alpha + gamma + d) + R*d;
    
ev(1,:) = [0, ev1/RT];
ev(2,:) = [ev(1,2), ev(1,2) + d*S/RT];
ev(3,:) = [ev(2,2), ev(2,2) + (beta_F*S*(I+epsilon*E)/N + beta_D*S*(I+epsilon*E))/RT];
ev(4,:) = [ev(3,2), ev(3,2) + d*E/RT];
ev(5,:) = [ev(4,2), ev(4,2) + kappa*E/RT];
ev(6,:) = [ev(5,2), ev(5,2) + I*(alpha+d)/RT];
ev(7,:) = [ev(6,2), ev(6,2) + gamma*I/RT];
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