%%% Stochastic iterations of the SIS model with 
%%% different set values of disease induced mortality alpha and baseline
%%% other parameters

%%% Used to collect data for Figures:
%%%     Fig 4.1

%%% Modify alpha for correct value of disease induced mortality

%%% IMPORTANT: this script uses the function STOCSIS. Inputs within this
%%% function need adapting when varying parameters

%%% would also recommend the use of a server. With a personal computer
%%% these iterations can take a very long time.

clearvars

K=1000;

T = 50; %end time
alpha = 10;

y0 = K*[0.995, 0.005];  %initial densities
sus = zeros(100,1); 
inf = zeros(100,1);
tvec = zeros(100,1);
sus(:,1) = y0(1);
inf(:,1) = y0(2);
tvec(:,1) = 0;

for i = 1:100

    t = 0;  %initial time
    y = y0;
    count = 2;
    
    while t < T
        res = STOCSIS(t, y, alpha);
        t = res(1);
        y = [res(2),res(3)];
        tvec(i, count) = t;
        sus(i, count) =  res(2);
        inf(i, count) = res(3);
        count = count + 1;
    end
end

save('StochAlpha10')