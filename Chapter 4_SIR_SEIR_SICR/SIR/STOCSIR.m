%%% Stochastic SIR model using Gillespie's Algorithm

function pop = STOCSIR(t, y, alpha, beta_D)
N = sum(y);
S = y(1);
I = y(2);
R = y(3);


%% Parameters
b = 10;
d = 1;
K = 1000;
q = (b - d)/(b*K);
%beta_D = 0.15;
beta_F = 30;

gamma = 50;
%alpha = 50;

%% MC Run

ev = zeros(6,2);  %event vector
tv = zeros(6,3);  %transition vector

ev1 = max(0, b*N*(1-q*N));
RT = ev1 + beta_F*S*I/N + beta_D*S*I + d*S + I*(alpha + gamma + d) + R*d;
    
ev(1,:) = [0, ev1/RT];
ev(2,:) = [ev(1,2), ev(1,2) + d*S/RT];
ev(3,:) = [ev(2,2), ev(2,2) + (beta_F*S*I/N + beta_D*S*I)/RT];
ev(4,:) = [ev(3,2), ev(3,2) + I*(alpha+d)/RT];
ev(5,:) = [ev(4,2), ev(4,2) + gamma*I/RT];
ev(6,:) = [ev(5,2), 1];

%set the upper limit to be greater than 1 so that the code is more
%efficient later on in choosing an event. It doens't matter what value this
%takes as all uniform random numbers lie between 0 and 1 inclusive.

tv(1,:) = [ 1, 0, 0];
tv(2,:) = [-1, 0, 0];
tv(3,:) = [-1, 1, 0];
tv(4,:) = [ 0,-1, 0];
tv(5,:) = [ 0,-1, 1];
tv(6,:) = [ 0, 0,-1];
 
U = rand(1);
for i = 1:6
    if U >= ev(i,1) && U < ev(i,2)
        y = y + tv(i,:);
        break
    end
end

dt = exprnd(1/RT);
t = t + dt;
 
pop = [t, y];

end