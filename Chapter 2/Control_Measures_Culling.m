%  Chapter 2: Varying Culling Control Measures Figures:
%       Fig 2.2, A.7, A.11, A.15
%    Running model simulations with control measures for:
%       Estonia under natural conditions Fig 2.2,
%       Estonia with supplemented feeding Fig A.7,
%       Spain under natural conditions Fig A.11,
%       Spain with supplemented feeding Fig A.15.

clearvars
options = odeset('Refine',1, 'RelTol',1e-9,'AbsTol',1e-9, 'NonNegative',1:8);
%Ode solver tolerances and restrictionsw

Scen = 1;       %distinguish between the scenarios
                       %Est natural conditions:         Scen = 1
                       %Est supplemented feeding:       Scen = 2
                       %Spain natural conditions:       Scen = 2.5
                       %Spain supplemented feeding:     Scen = 5
%Time of First Discovery tD, with outbreak starting in summer (t=0.6)
    %Estonia No Feeding tD = 0.7390
    %Estonia Feeding tD = 0.7261
    %Spain No Feeding tD = 1.4086
    %Spain Feeding tD = 0.6355
if Scen == 1
    implement = 0.1390;
elseif Scen == 2
    implement = 0.1261;
elseif Scen == 2.5
    implement = 0.8086;
elseif Scen == 5
    implement = 0.0355;
end

figure
for r=0:26:52
    TotDens = [];
    InfChrn = [];
    Carcass = [];
    for bC = 0:25
        b_C = bC/10;
        y0 = Scen*[0.796, 1.2, 0.002, 0.002, 0, 0, 0, 0]; %Initial Conditions
        %choose the relevant model formulation with or without feeding
        [t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, 0, 0), [0.6 0.6+implement], y0, options);
        %[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, 0, 0), [0.6 0.6+implement], y0, options);
        %[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, 0, 0), [0.6 0.6+implement], y0, options);
        %[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, 0, 0), [0.6 0.6+implement], y0, options);
        y0 = y(end,:);
        %choose the relevant model formulation with or without feeding
        [t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, b_C, r), [0.6+implement  5.6], y0, options);
        %[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, b_C, r), [0.6+implement  5.6], y0, options);
        %[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, b_C, r), [0.6+implement  5.6], y0, options);
        %[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, b_C, r), [0.6+implement  5.6], y0, options);
        totaldensity = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
        infected = y(:,3) + y(:,4);
        chronic = y(:,5) + y(:,6);
        
        f = find(t >= 2.6);
        g = find(t >= 3.6);
        totdens = mean(totaldensity(f(1):g(1)));
        infchrn = mean(infected(f(1):g(1))) + mean(chronic(f(1):g(1)));
        carcass = mean(y(f(1):g(1),7));
        
        TotDens = [TotDens, totdens];
        InfChrn = [InfChrn, infchrn];
        Carcass = [Carcass, carcass];
    end
    
    subplot(3,1,1)
    if r == 0
        plot(0:0.1:2.5, TotDens, 'k-');
    elseif r == 26
        plot(0:0.1:2.5, TotDens, 'k-.');
    else
        plot(0:0.1:2.5, TotDens, 'k:');
    end
    ylabel('Total Density')
    set(gca,'Xticklabel',[]);
    hold on

    subplot(3,1,2)
    if r == 0
        plot(0:0.1:2.5, InfChrn, 'k-');
    elseif r == 26
        plot(0:0.1:2.5, InfChrn, 'k-.');
    else
        plot(0:0.1:2.5, InfChrn, 'k:');
    end
    ylabel({'Inf. and'; 'Chron. Density'})
    set(gca,'Xticklabel',[]);
    hold on

    subplot(3,1,3)
    if r == 0
        plot(0:0.1:2.5, Carcass, 'k-');
    elseif r == 26
        plot(0:0.1:2.5, Carcass, 'k-.');
    else
        plot(0:0.1:2.5, Carcass, 'k:');
    end
    ylabel({'Infected';'Carcass Density'})
    hold on
end

subplot(3,1,1)
set(gca,'Position',[0.12 0.7 0.85 0.25])
str = 'A';
dim = [.87 .85 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,1,2)
%ylim([0 0.07])
set(gca,'Position',[0.12 0.4 0.85 0.25])
str = 'B';
dim = [.87 .55 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');
legend('r = 0', 'r = 26', 'r = 52', 'Location', 'North')

subplot(3,1,3)
% ylim([0 0.025])
xlabel('Value of b_C')
set(gca,'Position',[0.12 0.1 0.85 0.25])
str = 'C';
dim = [.87 .25 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');
