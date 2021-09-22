%%% Script to run the TB and ASF model simulations for Spain
%%% Can be used for the following figures in Chapter 3 of the thesis:
%%%     Fig 3.2 - 3.7, Fig B.1 and Fig B.2

%%% Uses the model function: ModelTBandASF

%%% Uses and explicit Euler Method to solve ODEs. This is done as we force
%%% the ASF infected population to 0 when the population density is very
%%% low. Changing population densities within a higher order solver is more
%%% problematic.

%%% 5 inputs are required here: 
%%%     P - prevalence level, 0 is low, 1 is high
%%%     d1 - carcass degradaion rate, 52/1 is high, 52/4 is low
%%%     tend - time to run simulations for
%%%     b_C - culling rate
%%%     r -  carcass removal rate

%%% time steps can be modified accordingly

clearvars

%% Model Inputs:

P = 0; 
%P=0 represents the low prevalence regions
%P=1 represents the high prevalence regions
tend = 26;
%tend denotes how long we run model for: tend - 1 years
%Control:
b_C = 0;   %Culling            v1: b_C = 0.2, r = 200
r = 0;     %Carcass Removal     v2: b_C = 0.7, r = 26
%High or Low Degradation Rate:
d1 = 52/4;  
%d1 = 52/1 represents high degradation rate
%d1 = 52/4 represents low degradation rate

%% Running the Model and Results

dt = 0.001;
t0 = 0;
tf = 1;
n = 1;

%%% Old Steady States: %y(:,1) = [1.2056, 0.0734, 0.0228, 0, 0, 0, 0, 0, 0, 0.7763, 0.0754, 0.0299, 0, 0, 0, 0, 0, 0, 1.6186, 0.1215, 0.0769, 0, 0, 0, 0, 0, 0, 0,0.0216];
if P == 0
    %Low Prevalence TB Steady State 
    y(:,1) = [1.2090, 0.0756, 0.0235, 0, 0, 0, 0, 0, 0, 0.7770, 0.0775, 0.0308, 0, 0, 0, 0, 0, 0, 1.6169, 0.1249, 0.0790, 0, 0, 0, 0, 0, 0, 0, 0.0222];
elseif P == 1
    %High Prevalence TB Steady State
    y(:,1) = [2.7877, 1.3761, 1.2843, 0, 0, 0, 0, 0, 0, 0.9310, 0.8192, 1.1964, 0, 0, 0, 0, 0, 0, 0.9698, 1.1904, 1.3792, 0, 0, 0, 0, 0, 0, 0, 0.6433];
end

%%% Old Steady States: %y(:,1) = [2.7750, 1.3911, 1.2983, 0, 0, 0, 0, 0, 0, 0.9207, 0.8250, 1.2069, 0, 0, 0, 0, 0, 0, 0.9516, 1.1935, 1.3879, 0, 0, 0, 0, 0, 0, 0,0.6488];


%Euler's method (1 year of endemic TB no ASF)
for t = t0+dt:dt:tf
    y(:,n+1) = y(:,n) + dt*TBandASFDDBwControl(t, y(:,n), P, d1, 0, 0);
    n = n + 1;
    for i = 1:29
        if y(i,end) <= 0
            y(i,end) = 0;
        end
    end
end

if P == 0
    %Low Prevalence introduction of ASF
    y(:,end+1) = [1.1970, 0.0749, 0.0233, 0.0120, 0.0007, 0.0002, 0, 0, 0, 0.7693, 0.0767, 0.0305, 0.0077, 0.0008, 0.0003, 0, 0, 0, 1.6009, 0.1237, 0.0782, 0.0160, 0.0012, 0.0008, 0, 0, 0, 0,0.0222];
  elseif P == 1
    %High Prevalence introduction of ASF
    y(:,end+1) = [2.7601, 1.3625, 1.2716, 0.0276, 0.0136, 0.0127, 0, 0, 0, 0.9218, 0.8111, 1.1846, 0.0092, 0.0081, 0.0118, 0, 0, 0, 0.9602, 1.1786, 1.3655, 0.0096, 0.0118, 0.0137, 0, 0, 0, 0, 0.6433];
 end
       
t0 = 1 + dt;
tf = tend;
n = n + 1;

for t = t0:dt:tf
    %%%%Euler's Method
    y(:,n+1) = y(:,n) + dt*TBandASFDDBwControl(t, y(:,n), P, d1, b_C, r);
    for i = 1:29
        if y(i,end) <= 0
            y(i,end) = 0;
        end
    end
    n = n + 1;
    
    %%%% Populations
    SASF = y(1,end)+y(2,end)+y(3,end)+y(10,end)+y(11,end)+y(12,end)+y(19,end)+y(20,end)+y(21,end);
    IASF = y(4,end)+y(5,end)+y(6,end)+y(13,end)+y(14,end)+y(15,end)+y(22,end)+y(23,end)+y(24,end);
    CASF = y(7,end)+y(8,end)+y(9,end)+y(16,end)+y(17,end)+y(18,end)+y(25,end)+y(26,end)+y(27,end);
    
    STB = y(1,end)+y(4,end)+y(7,end)+y(10,end)+y(13,end)+y(16,end)+y(19,end)+y(22,end)+y(25,end);
    
    N = SASF + IASF + CASF;
   
    if P == 0    
%%%%    %Low Prevalence%    %%%%
        if N - SASF <= 0.004  %%%% Artifically kill ASF if it is too low
            y(4,end)=0;y(7,end)=0;y(13,end)=0;y(16,end)=0;y(22,end)=0;y(25,end)=0;
            y(5,end)=0;y(8,end)=0;y(14,end)=0;y(17,end)=0;y(23,end)=0;y(26,end)=0;
            y(6,end)=0;y(9,end)=0;y(15,end)=0;y(18,end)=0;y(24,end)=0;y(27,end)=0;
            y(28,end)=0;
            b_C = 0;
            r = 0;
        end
        if N - STB <= 0.004  %%%%Artificially kill TB if it is too low
            y(2,end)=0;y(6,end)=0;y(11,end)=0;y(15,end)=0;y(20,end)=0;y(24,end)=0;
            y(3,end)=0;y(8,end)=0;y(12,end)=0;y(17,end)=0;y(21,end)=0;y(26,end)=0;
            y(5,end)=0;y(9,end)=0;y(14,end)=0;y(18,end)=0;y(23,end)=0;y(27,end)=0;
            y(29,end)=0;
        end
        if N <= 0.004  %%%%Artificially kill the population if it is too low
            y(1,end)=0;y(5,end)=0;y(9,end)=0;y(13,end)=0;y(17,end)=0;y(21,end)=0;
            y(2,end)=0;y(6,end)=0;y(10,end)=0;y(14,end)=0;y(18,end)=0;y(22,end)=0;
            y(3,end)=0;y(7,end)=0;y(11,end)=0;y(15,end)=0;y(19,end)=0;y(23,end)=0;
            y(4,end)=0;y(8,end)=0;y(12,end)=0;y(16,end)=0;y(20,end)=0;y(24,end)=0;
            y(25,end)=0;y(26,end)=0;y(27,end)=0;y(28,end)=0;y(29,end)=0;
        end
    elseif P == 1       
%%%%    %High Prevalence%    %%%%
        if N - SASF <= 0.012  
            y(4,end)=0;y(7,end)=0;y(13,end)=0;y(16,end)=0;y(22,end)=0;y(25,end)=0;
            y(5,end)=0;y(8,end)=0;y(14,end)=0;y(17,end)=0;y(23,end)=0;y(26,end)=0;
            y(6,end)=0;y(9,end)=0;y(15,end)=0;y(18,end)=0;y(24,end)=0;y(27,end)=0;
            y(28,end)=0;
            b_C = 0;
            r = 0;
        end
        if N - STB <= 0.012
            y(2,end)=0;y(6,end)=0;y(11,end)=0;y(15,end)=0;y(20,end)=0;y(24,end)=0;
            y(3,end)=0;y(8,end)=0;y(12,end)=0;y(17,end)=0;y(21,end)=0;y(26,end)=0;
            y(5,end)=0;y(9,end)=0;y(14,end)=0;y(18,end)=0;y(23,end)=0;y(27,end)=0;
            y(29,end)=0;
        end
        if N <= 0.012
            y(1,end)=0;y(5,end)=0;y(9,end)=0;y(13,end)=0;y(17,end)=0;y(21,end)=0;
            y(2,end)=0;y(6,end)=0;y(10,end)=0;y(14,end)=0;y(18,end)=0;y(22,end)=0;
            y(3,end)=0;y(7,end)=0;y(11,end)=0;y(15,end)=0;y(19,end)=0;y(23,end)=0;
            y(4,end)=0;y(8,end)=0;y(12,end)=0;y(16,end)=0;y(20,end)=0;y(24,end)=0;
            y(25,end)=0;y(26,end)=0;y(27,end)=0;y(28,end)=0;y(29,end)=0;
        end
    end
end

%% The Populations

totalpig = y(1,:)+y(2,:)+y(3,:)+y(4,:)+y(5,:)+y(6,:)+y(7,:)+y(8,:)+y(9,:);
totalyea = y(10,:)+y(11,:)+y(12,:)+y(13,:)+y(14,:)+y(15,:)+y(16,:)+y(17,:)+y(18,:);
totaladu = y(19,:)+y(20,:)+y(21,:)+y(22,:)+y(23,:)+y(24,:)+y(25,:)+y(26,:)+y(27,:);
totalpop = totalpig + totalyea + totaladu;

totalsusTB = y(1,:)+y(4,:)+y(7,:)+y(10,:)+y(13,:)+y(16,:)+y(19,:)+y(22,:)+y(25,:);
totalinfTB = y(2,:)+y(5,:)+y(8,:)+y(11,:)+y(14,:)+y(17,:)+y(20,:)+y(23,:)+y(26,:);
totalgenTB = y(3,:)+y(6,:)+y(9,:)+y(12,:)+y(15,:)+y(18,:)+y(21,:)+y(24,:)+y(27,:);

totalsusASF = y(1,:)+y(2,:)+y(3,:)+y(10,:)+y(11,:)+y(12,:)+y(19,:)+y(20,:)+y(21,:);
totalinfASF = y(4,:)+y(5,:)+y(6,:)+y(13,:)+y(14,:)+y(15,:)+y(22,:)+y(23,:)+y(24,:);
totalchrASF = y(7,:)+y(8,:)+y(9,:)+y(16,:)+y(17,:)+y(18,:)+y(25,:)+y(26,:)+y(27,:);

%% The Plots
tv = 0:dt:tf+dt;

%TB Densities Plots
figure
subplot(3,1,1)
plot(tv-1, totalpop, 'b-', tv-1, totalsusTB, 'g-', tv-1, totalinfTB, 'm-', tv-1, totalgenTB, 'r-', tv-1, totalinfTB+totalgenTB, 'k-')
ylabel('TB Dens.')
if P==0
    ylim([0 4.5])    %ylimits for low prevalence case
elseif P==1
    ylim([0 16])    %ylimits for high prevalence case
end
xlim([-1 tend-1])
set(gca, 'Xticklabels',[])
xticks([0 1 2 3 4 5 10 15 20 25])
%legend('Total','TB Susceptible','TB Infected','TB Generalised','TB Inf.+Gen.');
set(gca,'Position',[0.08 0.7 0.9 0.27])
str = 'D';
dim = [.9 .86 .03 .03];
annotation('textbox', dim, 'String', str, 'EdgeColor', 'w');

%ASF Densities plots
subplot(3,1,2)
plot(tv-1, totalpop, 'b-', tv-1, totalsusASF, 'g-', tv-1, totalinfASF, 'm-', tv-1, totalchrASF, 'r-', tv-1, totalinfASF+totalchrASF, 'k-')
ylabel('ASF Dens.')
if P==0
    ylim([0 4.5])    %ylimits for low prevalence case
elseif P==1
    ylim([0 16])    %ylimits for high prevalence case
end
xlim([-1 tend-1])
xticks([0 1 2 3 4 5 10 15 20 25])
%legend('Total','ASF Susceptible','ASF Infected','ASF Chronic','ASF Inf.+Chr.');
set(gca, 'Xticklabels',[])
set(gca,'Position',[0.08 0.4 0.9 0.27])
str = 'E';
dim = [.9 .55 .03 .03];
annotation('textbox', dim, 'String', str, 'EdgeColor', 'w');

%Prevalence of TB
subplot(3,1,3)
plot(tv-1, 100*(totalinfTB+totalgenTB)./totalpop, 'k-', tv-1, 100*totalinfTB./totalpop, 'k:', tv-1, 100*totalgenTB./totalpop, 'k-.')
ylabel('TB Prev. %')
if P==0
    ylim([0 10.5])    %ylimits for low prevalence case
elseif P==1
    ylim([0 65])     %ylimits for high prevalence case
end
xlim([-1 tend-1])
xlabel('Time t')
xticks([0 1 2 3 4 5 10 15 20 25])
%legend('Total','Infected','Generalised');
set(gca,'Position',[0.08 0.1 0.9 0.27])
str = 'F';
dim = [.9 .28 .03 .03];
annotation('textbox', dim, 'String', str, 'EdgeColor', 'w');


% totalpop(end)
% min(100*(totalinfTB+totalgenTB)./totalpop)
% min(totalinfTB)
% min(totalgenTB)

    



