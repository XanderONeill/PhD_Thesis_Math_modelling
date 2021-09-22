%%% Stochastic runs calculating extinction times for the SIS model with 
%%% different set values of disease induced mortality alpha and varying
%%% other parameters \gamma, \beta_F and \beta_D

%%% Both population extinction and disease extinction are considered
%%% In all scenarios an outbreak, defined as the susceptible population
%%% dropping below 900, must occur.

%%% Used to collect data for Figures:
%%%     Fig 4.3, and 4.6

%%% Modify alpha for correct value of disease induced mortality
%%% Within for loop change the parameter vector you would like to loop
%%% through to obtain data sets for that parameter. Change the ode solver
%%% inputs and the save command accordingly

%%% IMPORTANT: this script uses the function STOCSIS. Inputs within this
%%% function need adapting when varying parameters

%%% would also recommend the use of a server. With a personal computer
%%% these iterations can take a very long time.

clearvars
K = 1000;               %carrying capacity
y0 = K*[0.995,0.005];   %initial conditions

alpha = 80;             %disease induced mortality rate

betaDvec = (0:6)/20;    %dens. dep. trans. rate
betaFvec = 0:10:60;     %freq. dep. trans. rate
gammavec = 5:5:80;      %infected recovery rate

iter = 100;             %number of stochastic iterations

minorextinct = zeros(iter,length(gammavec));  %extinction with outbreak 
T = 500;                %length of time to run

for j = 1:length(gammavec)  
    gamma = gammavec(j);
    parfor i = 1:iter
        t = 0;
        y = y0;
        minsus = y0(1);
        while y(2) > 0 
            if t > T    %only calculate time to extinction below T = 50 years
                break
            end
            res = STOCSIS(t, y, alpha, gamma);
            t = res(1);
            y = [res(2),res(3)];
            minsus = min(minsus,y(1));
        end
        if minsus < 900 %outbreak occurs
            minorextinct(i,j) = t;
        end
    end
end
toc

save('SISgammaAlpha10it100t500');
