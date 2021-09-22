%%% Stochastic runs calculating extinction times for the SICS model with 
%%% different set values of disease induced mortality alpha and varying
%%% other parameters \epsilon and c (when epsilon >0). 

%%% This script is basically an extenstion of
%%% SICS_datasets. However, when the chronic population is infectious there
%%% are additional requirements for the pathogen to be considered extinct.
%%% Hence, the extra if loop within the code.

%%% Both population extinction and disease extinction are considered
%%% In all scenarios an outbreak, defined as the susceptible population
%%% dropping below 900, must occur.

%%% Used to collect data for Figures:
%%%     Fig 4.5, and 4.6

%%% Modify alpha for correct value of disease induced mortality
%%% Within for loop change the parameter vector you would like to loop
%%% through to obtain data sets for that parameter. Change the function
%%% inputs and the save command accordingly

%%% IMPORTANT: this script uses the function STOCSICS. Inputs within this
%%% function need adapting when varying parameters

%%% would also recommend the use of a server. With a personal computer
%%% these iterations can take a very long time.

clearvars
K = 1000;
y0 = K*[0.995,0.005,0];

alpha = 40;     %disease induced mortality rate

epsilonvec = 0:0.05:1; %infectiousness of exposed
cvec = 0:0.05:1;    %ratio of additional mortality for chronic due to infection

iter = 100; %number of stochastic iterations
minorextinct = zeros(iter,length(epsilonvec));  %extinction otherwise with outbreak
T = 1000;

for j = 1:length(epsilonvec)
    epsilon = epsilonvec(j);
    parfor i = 1:iter
        t = 0;
        y = y0;
        minsus = y0(1);
        if epsilon == 0
            while y(2) > 0 
                if t > T    %only calculate time to extinction below T = 50 years
                    break
                end
                res = STOCSICS(t, y, alpha, epsilon);
                t = res(1);
                y = [res(2),res(3),res(4)];
                minsus = min(minsus,y(1));
            end
        elseif epsilon > 0
            while y(2) + y(3) > 0 
                if t > T    %only calculate time to extinction below T = 50 years
                    break
                end
                res = STOCSICS(t, y, alpha, epsilon);
                t = res(1);
                y = [res(2),res(3),res(4)];
                minsus = min(minsus,y(1));
            end
        end    
        if minsus < 900          %minor outbreak
            minorextinct(i,j) = t;
        end
    end
end

save('SICSc05epsilonAlpha40it100t1000');