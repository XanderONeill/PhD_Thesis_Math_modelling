%%% Script for the analysis section of Chapter 3 in thesis. Run for figures
%%% Fig 3.8 and Fig 3.9 in the main text.

%%% Two model functions are used here. 1) ModelTB, used to produce Fig 3.8
%%% where no ASF death occurs and 2) ModelTBwASFdeath, used to produce Fig
%%% 3.9, where ASF death does occur.

%%% Part 1: Low prevalence regions
%%% Part 2: High prevalence regions

%%% adjust x limits accordingly in plots to produce each figure

clearvars

options = odeset('Refine',1, 'RelTol',1e-6,'AbsTol',1e-6,'NonNegative',1:10);

%% Low Prevalence Regions

susvecL = [];
infvecL = [];
genvecL = [];
totalvecL = [];
Kvec = 0:0.1:6;
for K = Kvec
    y0 = K*[0.2009,0.0142,0.0044,0.1279,0.0145,0.0058,0.2636,0.02317,0.0147,0.0041];
    [t,y] = ode45(@(t,y)  ModelTB(t, y, 0, K), [0 1000], y0, options);
    %[t,y] = ode45(@(t,y)  ModelTBwASFdeath(t, y, 0, K), [0 1000], y0, options);
    sus = y(:,1)+y(:,4)+y(:,7);
    inf = y(:,2)+y(:,5)+y(:,8);
    gen = y(:,3)+y(:,6)+y(:,9);
    totalpop = sus + inf + gen;
    susvecL = [susvecL; sus(end)];
    infvecL = [infvecL; inf(end)];
    genvecL = [genvecL; gen(end)];
    totalvecL = [totalvecL; totalpop(end)];
end

%% High Prevalence Regions
susvecH = [];
infvecH = [];
genvecH = [];
totalvecH = [];
Kvec1 = 0:0.5:30;
for K = Kvec1
    y0 = K*[0.2009,0.0142,0.0044,0.1279,0.0145,0.0058,0.2636,0.02317,0.0147,0.0041];
    [t,y] = ode45(@(t,y)  ModelTB(t, y, 1, K), [0 1000], y0, options);
    %[t,y] = ode45(@(t,y)  ModelTBwASFdeath(t, y, 1, K), [0 1000], y0, options);
    sus = y(:,1)+y(:,4)+y(:,7);
    inf = y(:,2)+y(:,5)+y(:,8);
    gen = y(:,3)+y(:,6)+y(:,9);
    totalpop = sus + inf + gen;
    susvecH = [susvecH; sus(end)];
    infvecH = [infvecH; inf(end)];
    genvecH = [genvecH; gen(end)];
    totalvecH = [totalvecH; totalpop(end)];
end

%% Plotting the results
figure
subplot(1,2,1)
plot(totalvecL, infvecL, 'k-', totalvecL, genvecL, 'k-.')
legend('Inf.','Gen.', 'Location', 'NorthWest')
xlabel('Total Population Density')
xlim([0 3])
ylabel('Steady State Densities')
ylim([-0.01 0.3])
set(gca,'Position',[0.09 0.1 0.42 0.85])
str = 'A';
dim = [.44 .89 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(1,2,2)
plot(totalvecH, infvecH, 'k-', totalvecH, genvecH, 'k-.')
legend('Inf.','Gen.', 'Location', 'NorthWest')
xlabel('Total Population Density')
xlim([0 5])
set(gca,'Position',[0.56 0.1 0.42 0.85])
str = 'B';
dim = [.91 .89 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

