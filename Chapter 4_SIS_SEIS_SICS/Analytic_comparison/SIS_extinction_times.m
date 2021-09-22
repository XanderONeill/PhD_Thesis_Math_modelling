%%% Analytical calculation of the mean time to pathogen extinction of the
%%% SIS model discussed in Chapter 4 of the Thesis.

%%% Also plots the mean time to pathogen extinction for a varying disease
%%% induced mortality rate

%%% Used in conjuction with Chapter4_Figure_2_run for SIS, to
%%% produce Figure Fig C.1

%%% Run the Chapter4_Figure_2_run script followed by this

% Xander_SIS_model
% SIS model from Xander O'Neill paper
% Compute mean time to extinction from any initial state

tic

clear 
close all
global maxN

I0 = 5; % Initial number infected

b = 10; % Maximum birth rate
d = 1; % Natural death rate
K = 1000; % Carrying capacity
betaF = 30; % Frequency dependent transmission parameter
betaD = 0.15; % Density dependent transmission parameter
gamma = 50; % Recovery rate
%alpha = 80; % Disease-induced mortality rate

% R0 = (betaF + betaD*K)/(alpha+gamma+d)

q = (b-d)/(b*K);
maxN = 1000;%round(1/q);

%state=zeros(maxN+1,maxN+1);
% Set up labels for states, with i=0 states first
%counter=0;
%for i=0:maxN
%    for s=0:maxN-i
%        counter=counter+1;
%        state(s+1,i+1)=counter;
%    end
%end
%maxcounter=counter;

% Note above loop is not a very efficient way to do things.
% Better would be to do a bit of algebra to work out the function
% that maps (i,s) to "state", and then instead of two loops simply define
% state = [appropriate function]
% so as to give same result as the above loops, but quicker.
% As follows:

state=ones(maxN+1,1)*(1 + (maxN+(3/2))*(0:maxN) - (1/2)*(0:maxN).^2)+(0:maxN)'*ones(1,maxN+1);

% Note it's important that the states with i=0 are enumerated first, so
% that they can be easily deleted from the transition rate matrix Q later.

% Avoid the double loop by defining a function instead
% The above double loop maps (s,i) to array element state(s+1,i+1)
% The function statef instead maps (s,i) to statef(s,i)
%
% (i,s) = (0,0),...,(0,maxN) enumerate as 1,...,maxN+1
% (i,s) = (1,0),...,(1,maxN-1) enumerate as (maxN+1)+1,...(maxN+1)+maxN
% (i,s) = (2,0),...,(2,maxN-2) enumerate as (maxN+1)+maxN+1,....(maxN+1)+maxN+(maxN-1)
%
% (i,0) for general i enumerates as 1+(maxN+1)+(maxN)+(maxN-1)+...+(maxN+2-i)
% That is, 1 + \sum_{j=1}^i (maxN+2-j)
% That is, 1 + i*(maxN+2) - (1/2)*i*(i+1)
% That is, 1 + i*(maxN+(3/2)) - (1/2)*i^2
% So then
% (i,s) enumerates as 1 + (maxN+(3/2))*i - (1/2)*i^2 + s
% And
% maxcounter = 1 + (maxN+(3/2))*maxN - (1/2)*maxN^2 + 0 
% = 1 + (3/2)*maxN + (1/2)*maxN^2

maxcounter = 1 + (3/2)*maxN + (1/2)*maxN^2;

nonzero_rates = (5/2)*maxN^2 + (3/2)*maxN - 1; 
% See calculation below of how many nonzero off-diagonal entries there are
% in the transition rate matrix

alpha_set = 30:5:80;
mean_extinction_times = zeros(1,length(alpha_set));
R0_values = zeros(1,length(alpha_set));
outbreak_probs = zeros(1,length(alpha_set));

for alpha=alpha_set
    
rowindex=zeros(1,nonzero_rates);
columnindex=zeros(1,nonzero_rates);
entries=zeros(1,nonzero_rates);

counter=0;

s=0;
    for i=1:maxN-s-1
        counter=counter+1;
        rowindex(counter)=state(s+1,i+1);
        columnindex(counter)=state(s+1+1,i+1);
        entries(counter)=b*(s+i)*(1-q*(s+i));
        % Q(state(s+1,i+1),state(s+1+1,i+1)) = b*(s+i)*(1-q*(s+i)); % s=0:maxN-1, i=0:maxN-s-1
    end % Number of entries here is maxN-1

for s=1:maxN-1
    for i=0:maxN-s-1
        counter=counter+1;
        rowindex(counter)=state(s+1,i+1);
        columnindex(counter)=state(s+1+1,i+1);
        entries(counter)=b*(s+i)*(1-q*(s+i));
        % Q(state(s+1,i+1),state(s+1+1,i+1)) = b*(s+i)*(1-q*(s+i)); % s=0:maxN-1, i=0:maxN-s-1
    end 
end % Number of entries here is (1/2)*maxN*(maxN-1)

for s=1:maxN
    for i=0:maxN-s
        counter=counter+1;
        rowindex(counter)=state(s+1,i+1);
        columnindex(counter)=state(s-1+1,i+1);
        entries(counter)=d*s;
        % Q(state(s+1,i+1),state(s-1+1,i+1)) = d*s; % s=1:maxN
    end
end % Number of entries here is (1/2)*maxN*(maxN+1)

for s=1:maxN-1
    for i=1:maxN-s
        counter=counter+1;
        rowindex(counter)=state(s+1,i+1);
        columnindex(counter)=state(s-1+1,i+1+1);
        entries(counter)=betaF*s*i/(s+i) + betaD*s*i;
        % Q(state(s+1,i+1),state(s-1+1,i+1+1)) = betaF*s*i/(s+i) + betaD*s*i; % s=1:maxN-1, i=1:maxN-s  
    end
end % Number of entries here is (1/2)*maxN*(maxN-1)

for s=0:maxN-1
    for i=1:maxN-s
        counter=counter+1;
        rowindex(counter)=state(s+1,i+1);
        columnindex(counter)=state(s+1+1,i-1+1);
        entries(counter)=gamma*i;
        % Q(state(s+1,i+1),state(s+1+1,i-1+1)) = gamma*i; % s=0:maxN-1, i=1:maxN-s
    end
end % Number of entries here is (1/2)*maxN*(maxN+1) 

for s=0:maxN-1
    for i=1:maxN-s
        counter=counter+1;
        rowindex(counter)=state(s+1,i+1);
        columnindex(counter)=state(s+1,i-1+1);
        entries(counter)=(alpha+d)*i;
        % Q(state(s+1,i+1),state(s+1,i-1+1)) = (alpha+d)*i; % s=0:maxN-1, i=1:maxN-s
    end
end % Number of entries here is (1/2)*maxN*(maxN+1)

% Total number of nonzero entries is (5/2)*maxN^2 + (3/2)*maxN - 1

% last=length(nonzeros(columnindex)); % Check - evaluates to (5/2)*maxN^2+(3/2)*maxN-1

Q = sparse(rowindex(1:nonzero_rates),columnindex(1:nonzero_rates),entries(1:nonzero_rates),maxcounter,maxcounter,nonzero_rates);
Q = Q - diag(sum(Q')); % Fill in diagonal entries so that each row sums to zero

% Truncate by removing all i=0 states (corresponding to extinction)
QC = Q(maxN+2:end,maxN+2:end);

% [evect,eval]=eigs(QC',1,'LR'); % Compute eigenvalue corresponding to QSD
% QSD = evect/sum(evect); % Normalise QSD
% QSD_expected_extinction_time = -1/eval % Expected time to extinction starting from QSD

expected_extinction_times = -QC\ones(maxcounter-maxN-1,1); % Expected time to extinction starting from any initial state

% meantimes = zeros(1,K);
% for I0=1:K;
% meantimes(I0)=expected_extinction_times(state(K-I0+1,I0+1)-maxN-1);
% end
% figure
% plot(1:K,meantimes) % Times to extinction starting from population size S+I = K, as a function of I=1,...,K
% meantimes(5) % Time to extinction starting from I=5, S=K-5

mean_extinction_times(find(alpha_set == alpha)) = expected_extinction_times(state(K-I0+1,I0+1)-maxN-1); %#ok<FNDSB>
R0_values(find(alpha_set == alpha)) = (betaF + betaD*K)/(alpha+gamma+d);
outbreak_probs(find(alpha_set == alpha)) = 1-((alpha+gamma+d)/(betaF + betaD*K))^I0;

% Would be a good idea to save the whole of the expected_extinction_times
% vector, at the moment most of it is discarded.

end

R0_values
outbreak_probs

figure
plot(alpha_set,mean_extinction_times,'b.-')
hold on
xlim([20 80])
ylim([0 100])
ylabel('Mean time to pathogen extinction, \tau')
xlabel('Disease-induced mortality rate, \alpha')

toc