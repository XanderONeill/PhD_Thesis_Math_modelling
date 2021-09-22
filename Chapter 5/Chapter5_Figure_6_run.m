%%% Script to produce Figure 5.6 in the thesis, Chapter 5

%%% Uses the tick-hosts model with host specific density dependence,
%%% representative of equations (5.6) in the thesis

%%% Loops through the small mammal host density values and ungulate host
%%% density values to produce a 3D plot of tick densities

clearvars
options = odeset('Refine', 1,'RelTol',1e-5,'AbsTol',1e-5,'NonNegative',1:7);

Huvec = [1:5,10:5:100];
Hsvec = [1:10,20:5:200];
LQStar = zeros(length(Hsvec),length(Huvec)); LFStar = zeros(length(Hsvec),length(Huvec));
NQStar = zeros(length(Hsvec),length(Huvec)); NFStarS = zeros(length(Hsvec),length(Huvec));
                                             NFStarU = zeros(length(Hsvec),length(Huvec));
AQStar = zeros(length(Hsvec),length(Huvec)); AFStar = zeros(length(Hsvec),length(Huvec));

for i = 1:length(Hsvec)     %looping through the host densities
    for j = 1:length(Huvec)
        Hs = Hsvec(i);
        Hu = Huvec(j);
        
        y0 = 100*ones(7,1);
        [t,y] = ode15s(@(t,y) Model_tick_host(t, y, Hs, Hu), [0 10000], y0, options);
        
        LQStar(i,j) = y(end,1);    LFStar(i,j) = y(end,2);
        NQStar(i,j) = y(end,3);    NFStarS(i,j) = y(end,4);
                                   NFStarU(i,j) = y(end,5);
        AQStar(i,j) = y(end,6);    AFStar(i,j) = y(end,7);
    end
end

%% Figure 5.6: Class Questing and Feeding Steady states
figure
subplot(2,3,1)
surf(Hsvec, Huvec, LQStar')
set(gca, 'Position', [0.06 0.70 0.24 0.2520])
str = 'A(i)';
dim = [.24 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(2,3,2)
surf(Hsvec, Huvec, NQStar')
set(gca, 'Position', [0.38 0.70 0.24 0.2520])
str = 'B(i)';
dim = [.56 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');

subplot(2,3,3)
surf(Hsvec, Huvec, AQStar')
set(gca, 'Position', [0.7 0.70 0.24 0.2520])
str = 'C(i)';
dim = [.88 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');

subplot(2,3,4)
surf(Hsvec, Huvec, LFStar')
set(gca, 'Position', [0.06 0.4 0.24 0.2520])
str = 'A(ii)';
dim = [.24 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(2,3,5)
surf(Hsvec, Huvec, (NFStarS+NFStarU)')
set(gca, 'Position', [0.38 0.4 0.24 0.2520])
str = 'B(ii)';
dim = [.56 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');

subplot(2,3,6) 
surf(Hsvec, Huvec, AFStar')
set(gca, 'Position', [0.7 0.4 0.24 0.2520])
str = 'C(ii)';
dim = [.88 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');

