%%% Producing Figures 5.8, 5.9 and 5.10 in the thesis, Chapter 5

%%% Uses Model_tick_host_wInf
%%% Model simulations for the tick-host model with specific host density
%%% depedence and infection

%%% Produces 3 figures, one for the prevalence of ticks and hosts, one for
%%% the density of infected ticks and one for the density of each class of
%%% host (susceptible, infected and recovered)

clearvars
%close all

options = odeset('Refine', 1,'RelTol',1e-5,'AbsTol',1e-5,'NonNegative',1:20);
       %%% ode tolerances and restrictions
       
%% Running the model simulations with infection
Kuvec = [1:5,10:5:100];
Ksvec = [1:10,20:5:200];

LqS = zeros(length(Ksvec),length(Kuvec));
LqI = LqS; LfS = LqS; LfI = LqS; 
NqS = LqS; NqI = LqS; NfsS = LqS; NfsI = LqS; NfuS = LqS; NfuI = LqS;
AqS = LqS; AqI = LqS; AfS = LqS; AfI = LqS;

HsS = LqS; HsI = LqS; HsR = LqS;
HuS = LqS; HuI = LqS; HuR = LqS;

for i = 1:length(Ksvec)
    for j = 1:length(Kuvec)
        K_s = Ksvec(i);
        K_u = Kuvec(j);
        
        y0 = K_s*K_u*10*ones(20,1);
        y0(15:20) = [K_s,0,0,K_u,0,0];
        [t,y] = ode45(@(t,y) Model_tick_host_wInf(t, y, K_s, K_u), [0 1500000], y0, options);
        
        LqS(i,j) = y(end,1); LqI(i,j) = y(end,2); LfS(i,j) = y(end,3);
        LfI(i,j) = y(end,4); 
        NqS(i,j) = y(end,5); NqI(i,j) = y(end,6); NfsS(i,j) = y(end,7);
        NfsI(i,j) = y(end,8); NfuS(i,j) = y(end,9); NfuI(i,j) = y(end,10); 
        AqS(i,j) = y(end,11); AqI(i,j) = y(end,12); AfS(i,j) = y(end,13); 
        AfI(i,j) = y(end,14); 
        
        HsS(i,j) = y(end,15); HsI(i,j) = y(end,16); HsR(i,j) = y(end,17); 
        HuS(i,j) = y(end,18); HuI(i,j) = y(end,19); HuR(i,j) = y(end,20); 
    end
end

%% Determining the relevant population densities
totalticks = LqS + LqI + LfS + LfI + NqS + NqI +NfsS + NfsI + NfuS + NfuI + AqS + AqI + AfS + AfI;
infticks = LqI + LfI + NqI + NfsI + NfuI + AqI + AfI;

totalungulates = HuS + HuI + HuR;
seroungulates  = HuI + HuR;
totalsmall = HsS + HsI + HsR;
serosmall  = HsI + HsR;

infnym = NqI + NfsI + NfuI;
infadu = AqI + AfI;

TL = LqS+LqI+LfS+LfI;
TN = NqS+NqI+NfsS+NfuS+NfsI+NfuI;
TA = AqS+AqI+AfS+AfI;
        

%% Figure 5.8
%%% 4x3 Figure plotting prevalences with varying small mammal density and
%%% varying ungulate density
%%% Row 1: Questing larvae, nymph, adult
%%% Row 2: Feeding larvae, nymph, adult
%%% Row 3: Total larvae, nymph, adult
%%% Row 4: Total ticks, small mammals, ungulates

figure
subplot(4,3,1)
surf(Ksvec, Kuvec, 100*(LqI./(LqS+LqI))')
set(gca, 'Position', [0.04 0.79 0.18 0.18])
title('A(i): Ques. Lar.')
subplot(4,3,2)
surf(Ksvec, Kuvec, 100*(NqI./(NqS+NqI))')
set(gca, 'Position', [0.26 0.79 0.18 0.18])
title('A(ii): Ques. Nym.')
subplot(4,3,3)
surf(Ksvec, Kuvec, 100*(AqI./(AqS+AqI))')
set(gca, 'Position', [0.48 0.79 0.18 0.18])
title('A(iii): Ques. Adu.')
subplot(4,3,4)
surf(Ksvec, Kuvec, 100*(LfI./(LfS+LfI))')
set(gca, 'Position', [0.04 0.54 0.18 0.18])
title('B(i): Feed. Lar.')
subplot(4,3,5)
surf(Ksvec, Kuvec, 100*((NfsI+NfuI)./(NfsS+NfuS+NfsI+NfuI))')
set(gca, 'Position', [0.26 0.54 0.18 0.18])
title('B(ii): Feed. Nym.')
subplot(4,3,6)
surf(Ksvec, Kuvec, 100*(AfI./(AfS+AfI))')
set(gca, 'Position', [0.48 0.54 0.18 0.18])
title('B(iii): Feed. Adu.')
subplot(4,3,7)
surf(Ksvec, Kuvec, 100*((LfI+LqI)./TL)')
set(gca, 'Position', [0.04 0.29 0.18 0.18])
title('C(i): Tot. Lar.')
subplot(4,3,8)
surf(Ksvec, Kuvec, 100*((NfsI+NfuI+NqI)./TN)')
set(gca, 'Position', [0.26 0.29 0.18 0.18])
title('C(ii): Tot. Nym.')
subplot(4,3,9)
surf(Ksvec, Kuvec, 100*((AfI+AqI)./TA)')
set(gca, 'Position', [0.48 0.29 0.18 0.18])
title('C(iii): Tot. Adu.')
subplot(4,3,10)
surf(Ksvec, Kuvec, (100*infticks./totalticks)')
set(gca, 'Position', [0.04 0.04 0.18 0.18])
title('D(i): Tot. Ticks')
subplot(4,3,11)
surf(Ksvec, Kuvec, 100*((HsI+HsR)./(HsS+HsI+HsR))')
set(gca, 'Position', [0.26 0.04 0.18 0.18])
title('D(ii): Small Mam.')
subplot(4,3,12)
surf(Ksvec, Kuvec, 100*((HuI+HuR)./(HuS+HuI+HuR))')
set(gca, 'Position', [0.48 0.04 0.18 0.18])
title('D(iii): Ungulate')

%% Figure 5.9
%%% 2x3 figure plotting infected densities for
%%% Questing larvae, nymph and adults
%%% Feeding larvae, nymph and adults

figure
subplot(2,3,1)
surf(Ksvec, Kuvec, LqI')
set(gca, 'Position', [0.06 0.72 0.24 0.24])
title('A(i): Inf. Ques. Larvae')
subplot(2,3,2)
surf(Ksvec, Kuvec, NqI')
set(gca, 'Position', [0.385 0.72 0.24 0.24])
title('B(i): Inf. Ques. Nymph')
subplot(2,3,3)
surf(Ksvec, Kuvec, AqI')
set(gca, 'Position', [0.71 0.72 0.24 0.24])
title('C(i): Inf. Ques. Adult')
subplot(2,3,4)
surf(Ksvec, Kuvec, LfI')
set(gca, 'Position', [0.06 0.4 0.24 0.24])
title('A(ii): Inf. Feed. Larvae')
subplot(2,3,5)
surf(Ksvec, Kuvec, (NfsI+NfuI)')
set(gca, 'Position', [0.385 0.4 0.24 0.24])
title('B(ii): Inf. Feed. Nymph')
subplot(2,3,6)
surf(Ksvec, Kuvec, AfI')
set(gca, 'Position', [0.71 0.4 0.24 0.24])
title('C(ii): Inf. Feed. Adult')

%% Figure 5.10
%%% 2x3 figure plotting infected densities for
%%% Susceptible, Infected and Recovered small mammal hosts
%%% Susceptible, Infected and Recovered ungulate hosts

figure
subplot(2,3,1)
surf(Ksvec, Kuvec, HsS')
set(gca, 'Position', [0.06 0.71 0.24 0.24])
title('A(i): H_{Ss} Dens.')
subplot(2,3,2)
surf(Ksvec, Kuvec, HsI')
set(gca, 'Position', [0.385 0.71 0.24 0.24])
title('B(i): H_{Si} Dens.')
subplot(2,3,3)
surf(Ksvec, Kuvec, HsR')
set(gca, 'Position', [0.71 0.71 0.24 0.24])
title('C(i): H_{Sr} Dens.')
subplot(2,3,4)
surf(Ksvec, Kuvec, HuS')
set(gca, 'Position', [0.06 0.39 0.24 0.24])
title('A(ii): H_{Us} Dens.')
subplot(2,3,5)
surf(Ksvec, Kuvec, HuI')
set(gca, 'Position', [0.385 0.39 0.24 0.24])
title('B(ii): H_{Ui} Dens.')
subplot(2,3,6)
surf(Ksvec, Kuvec, HuR')
set(gca, 'Position', [0.71 0.39 0.24 0.24])
title('C(ii): H_{Ur} Dens.')
