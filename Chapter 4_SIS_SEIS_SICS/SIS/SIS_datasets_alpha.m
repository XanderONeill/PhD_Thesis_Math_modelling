%%% Stochastic runs calculating extinction times for the SIS model with 
%%% a varying level of disease induced mortality alpha for fixed
%%% other parameters \gamma, \beta_F and \beta_D

%%% Both population extinction and disease extinction are considered
%%% In all scenarios an outbreak, defined as the susceptible population
%%% dropping below 900, must occur.

%%% Used to collect data for Figures:
%%%     Fig 4.2

%%% IMPORTANT: this script uses the function STOCSIS. Inputs within this
%%% function need adapting when varying parameters

%%% Change number of iterations and time to stop for more/less accurate
%%% results

%%% would also recommend the use of a server. With a personal computer
%%% these iterations can take a very long time.



clearvars
K = 1000;
y0 = K*[0.995,0.005];

alphavec = 5:5:80;     %disease induced mortality rate
iter = 200; %number of stochastic iterations

minorextinct = zeros(iter,length(alphavec));  %extinction otherwise with outbreak
T = 1000;

for j = 1:length(alphavec)
    alpha = alphavec(j);
    parfor i = 1:iter
        t = 0;
        y = y0;
        minsus = y0(1);
        while y(2) > 0 
            if t > T    %only calculate time to extinction below T = 50 years
                break
            end
            res = STOCSIS(t, y, alpha);
            t = res(1);
            y = [res(2),res(3)];
            minsus = min(minsus,y(1));
        end    
        if minsus < 900   %total pop drops below 900 minor outbreak
            minorextinct(i,j) = t;  
        end
    end
end

save('SISAlphait100t1000');