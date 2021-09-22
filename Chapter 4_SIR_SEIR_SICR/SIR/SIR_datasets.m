%%% Stochastic runs calculating extinction times for the SIR model with 
%%% different set values of disease induced mortality alpha and varying
%%% other parameters \gamma, \beta_F and \beta_D

%%% Both population extinction and disease extinction are considered
%%% In all scenarios an outbreak, defined as the susceptible population
%%% dropping below 900, must occur.

%%% Used to collect data for Figures:
%%%     Fig 4.4, and 4.7

%%% Modify alpha for correct value of disease induced mortality
%%% Within for loop change the parameter vector you would like to loop
%%% through to obtain data sets for that parameter. Change the function
%%% inputs and the save command accordingly

%%% IMPORTANT: this script uses the function STOCSIR. Inputs within this
%%% function need adapting when varying parameters

%%% would also recommend the use of a server. With a personal computer
%%% these iterations can take a very long time.

clearvars
K = 1000;
y0 = K*[0.995,0.005,0];

alpha = 1;     %disease induced mortality rate

betaDvec = (0:6)/20; %dens. dep. trans. rate
betaFvec = 0:10:60; %freq. dep. trans. rate
gammavec = 5:5:80; %infected recovery rate

iter = 100; %number of stochastic iterations
T = 500;
 
minorextinct = zeros(iter,length(betaDvec));  %extinction otherwise with outbreak

for j = 1:length(betaDvec)
    betaD = betaDvec(j);
    parfor i = 1:iter
        t = 0;
        y = y0;
        minsus = y0(1);
        while y(2) > 0 
            if t > T    %only calculate time to extinction below T = 50 years
                break
            end
            res = STOCSIR(t, y, alpha, betaD);
            t = res(1);
            y = [res(2),res(3),res(4)];
            minsus = min(minsus,y(1));
        end
        if minsus < 900       %minor outbreak
            minorextinct(i,j) = t;
        end
    end
end

save('SIRbetaDAlpha1it100t500');