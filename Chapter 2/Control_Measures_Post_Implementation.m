%  Chapter 2: Varying Post Discovery Implementation Time
%  of combined control measures, where b_C = 0.75 and r = 52. 
%  Representative Figures:
%       Fig 2.4, A.9, A.13, A.17
%    Running model simulations with control measures for:
%       Estonia under natural conditions Fig 2.4,
%       Estonia with supplemented feeding Fig A.9,
%       Spain under natural conditions Fig A.13,
%       Spain with supplemented feeding Fig A.17.

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
for event = 1:2
    TotDens = [];
    InfChrn = [];
    Carcass = [];
    if event == 1               %implement control immediately
        implement = 0.1390;
    else                        %delay implementation of control by 6 months
        implement = 0.7390;
    end
    y0 = Scen*[0.796, 1.2, 0.002, 0.002, 0, 0, 0, 0]; %Initial Conditions
    yfini = [];
    tfini = [];
    %choose the relevant model formulation with or without feeding
    [t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, 0, 0), [0.6 0.6+implement], y0, options);
    %[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, 0, 0), [0.6 0.6+implement], y0, options);
    %[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, 0, 0), [0.6 0.6+implement], y0, options);
    %[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, 0, 0), [0.6 0.6+implement], y0, options);
    y0 = y(end,:);
    yfini = [yfini; y];
    tfini = [tfini; t];
    %choose the relevant model formulation with or without feeding
    [t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, 0.75, 52), [0.6+implement  5.6], y0, options);
    %[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, 0.75, 52), [0.6+implement  5.6], y0, options);
    %[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, 0.75, 52), [0.6+implement  5.6], y0, options);
    %[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, 0.75, 52), [0.6+implement  5.6], y0, options);
    yfini = [yfini; y];
    tfini = [tfini; t];
    
    totdens = yfini(:,1) + yfini(:,2) + yfini(:,3) + yfini(:,4) + yfini(:,5) + yfini(:,6);
    infected = yfini(:,3) + yfini(:,4);
    chronic = yfini(:,5) + yfini(:,6);
        
    subplot(3,1,1)
    if event == 1
        plot(tfini - 0.6, totdens, 'k-');
    else
        plot(tfini - 0.6, totdens, 'k-.');
    end
    set(gca,'Xticklabel',[]);
    ylabel('Total Density')
    hold on

    subplot(3,1,2)
    if event == 1
        plot(tfini - 0.6, infected + chronic, 'k-');
    else
        plot(tfini - 0.6, infected + chronic, 'k-.');
    end
    set(gca,'Xticklabel',[]);
    ylabel({'Inf. and'; 'Chron. Density'})
    hold on

    subplot(3,1,3)
    if event == 1
        plot(tfini - 0.6, yfini(:,7), 'k-');
    else
        plot(tfini - 0.6, yfini(:,7), 'k-.');
    end
    ylabel({'Infected';'Carcass Density'})
    hold on
end

subplot(3,1,1)
%ylim([0 4.5])
yticks([0 2.5 5 7.5 10])
set(gca,'Position',[0.12 0.7 0.85 0.25])
str = 'A';
dim = [.87 .85 .05 .08];
annotation('textbox', dim, 'String', str, 'EdgeColor', 'w');

subplot(3,1,2)
%ylim([0 0.41])
set(gca,'Position',[0.12 0.4 0.85 0.25])
str = 'B';
dim = [.87 .55 .05 .08];
annotation('textbox', dim, 'String', str, 'EdgeColor', 'w');
legend('Discovered','6 months later', 'Location', 'North')

subplot(3,1,3)
% ylim([0 0.5])
% yticks([0 0.1 0.2 0.3 0.4 0.5])
xlabel('Time')
set(gca,'Position',[0.12 0.1 0.85 0.25])
str = 'C';
dim = [.87 .25 .05 .08];
annotation('textbox', dim, 'String', str, 'EdgeColor', 'w');