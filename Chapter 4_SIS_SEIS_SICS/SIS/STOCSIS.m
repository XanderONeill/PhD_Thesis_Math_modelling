function pop = STOCSIS(t, y, alpha)
N = sum(y);
S = y(1);
I = y(2);

%% Parameters
b = 10;
d = 1;
K = 1000;
q = (b - d)/(b*K);
beta_D = 0.15;  %0.15;
beta_F = 30; %beta_D*K/5;

gamma = 50;
%alpha =50;

%% MC Run
ev = zeros(5,2);  %event vector
tv = zeros(5,2);  %transition vector

ev1 = max(0, b*N*(1-q*N));
RT = ev1 + beta_F*S*I/N + beta_D*S*I + d*S + gamma*I + I*(alpha + d);
 
ev(1,:) = [0, ev1/RT];
ev(2,:) = [ev(1,2), ev(1,2) + (beta_F*S*I/N + beta_D*S*I)/RT];
ev(3,:) = [ev(2,2), ev(2,2) + d*S/RT];
ev(4,:) = [ev(3,2), ev(3,2) + gamma*I/RT];
ev(5,:) = [ev(4,2), 1]; 

tv(1,:) = [ 1, 0];
tv(2,:) = [-1, 1];
tv(3,:) = [-1, 0];
tv(4,:) = [ 1,-1];
tv(5,:) = [ 0,-1];
 
U = rand(1);
for i = 1:5
    if U >= ev(i,1) && U < ev(i,2)
        y = y + tv(i,:);
        break
    end
end

dt = exprnd(1/RT);
t = t + dt;
 
pop = [t, y];

end