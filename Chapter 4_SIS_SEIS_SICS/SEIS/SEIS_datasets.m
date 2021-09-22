%%% Stochastic runs calculating extinction times for the SIS model with 
%%% different set values of disease induced mortality alpha and varying
%%% other parameters \gamma, \beta_F, \beta_D, \kappa and \epsilon

%%% Both population extinction and disease extinction are considered
%%% In all scenarios an outbreak, defined as the susceptible population
%%% dropping below 900, must occur.

%%% Used to collect data for Figures:
%%%     Fig 4.3

%%% Modify alpha for correct value of disease induced mortality
%%% Within for loop change the parameter vector you would like to loop
%%% through to obtain data sets for that parameter. Change the function
%%% inputs and the save command accordingly

%%% IMPORTANT: this script uses the function STOCSEIS. Inputs within this
%%% function need adapting when varying parameters

%%% would also recommend the use of a server. With a personal computer
%%% these iterations can take a very long time.

tic
clearvars
K = 1000;
y0 = K*[0.995,0,0.005];

alpha = 80;     %disease induced mortality rate

betaDvec = (0:6)/20; %dens. dep. trans. rate
betaFvec = 0:10:60; %freq. dep. trans. rate
gammavec = 5:5:80; %infected recovery rate
kappavec = 10:10:160; %exposed recovery rate
epsilonvec = 0:0.05:1; %infectiousness of exposed

iter = 100; %number of stochastic iterations

minorextinct = zeros(iter,length(kappavec));  %extinction otherwise with outbreak
T = 1000;

for j = 1:length(kappavec)
    kappa = kappavec(j);
    parfor i = 1:iter
        t = 0;
        y = y0;
        minsus = y0(1);
        while y(2) + y(3) > 0 
            if t > T    %only calculate time to extinction below T = 500 years
                break
            end
            res = STOCSEIS(t, y, alpha, kappa);
            t = res(1);
            y = [res(2),res(3),res(4)];
            minsus = min(minsus,y(1));
        end
        if minsus < 900       %outbreak
            minorextinct(i,j) = t;
        end
    end
end
toc

save('SEISkappaAlpha80it100t500');
