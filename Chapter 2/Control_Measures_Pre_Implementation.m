%  Chapter 2: Implementation Control Measures Figures:
%       Fig 2.5, A.10, A.14, A.18
%    Running model simulations with control measures for:
%       Estonia under natural conditions Fig 2.5,
%       Estonia with supplemented feeding Fig A.10,
%       Spain under natural conditions Fig A.14,
%       Spain with supplemented feeding Fig A.18.

clearvars

options = odeset('Refine',1, 'RelTol',1e-8,'AbsTol',1e-9, 'NonNegative',1:8);
%Ode solver tolerances and restrictions

b_C = 0.75; r = 52; %combination of control measures

Scen = 2;       %distinguish between the scenarios
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

%% Implementation 1 year before start of infection 

y0 = Scen*[0.798, 1.202, 0, 0, 0, 0, 0, 0]; %Initial Conditions

%choose the relevant model formulation with or without feeding
%[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, b_C, r), [0.6 1.6], y0, options);
[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, b_C, r), [0.6 1.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, b_C, r), [0.6 1.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, b_C, r), [0.6 1.6], y0, options);
y0 = y(end,:);
y0(1) = y0(1)-Scen*0.002;
y0(2) = y0(2)-Scen*0.002;
y0(3) = Scen*0.002;
y0(4) = Scen*0.002;
y1 = y;
t1 = t;
%choose the relevant model formulation with or without feeding
%[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, b_C, r), [1.6  6.6], y0, options);
[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, b_C, r), [1.6 6.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, b_C, r), [1.6 6.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, b_C, r), [1.6 6.6], y0, options);
y = [y1; y];
t = [t1; t];
totalpop = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
infandchron = y(:,3) + y(:,4) + y(:,5) + y(:,6);
carcass = y(:,7);


figure  
subplot(3,1,1)
plot(t-0.6, totalpop, 'k:')
hold on
subplot(3,1,2)
plot(t-0.6, infandchron, 'k:')
hold on
subplot(3,1,3)
plot(t-0.6, carcass, 'k:')
hold on


%% Implementation 6 months before start of infection

y0 = Scen*[0.798, 1.202, 0, 0, 0, 0, 0, 0]; %Initial Conditions

%choose the relevant model formulation with or without feeding
%[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, 0, 0), [0.6 1.1], y0, options);
[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, 0, 0), [0.6 1.1], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, 0, 0), [0.6 1.1], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, 0, 0), [0.6 1.1], y0, options);
y0 = y(end,:);
y1 = y;
t1 = t;
%choose the relevant model formulation with or without feeding
%[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, b_C, r), [1.1 1.6], y0, options);
[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, b_C, r), [1.1 1.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, b_C, r), [1.1 1.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, b_C, r), [1.1 1.6], y0, options);
y0 = y(end,:);
y0(1) = y0(1)-Scen*0.002;
y0(2) = y0(2)-Scen*0.002;
y0(3) = Scen*0.002;
y0(4) = Scen*0.002;
y2 = y;
t2 = t;
%choose the relevant model formulation with or without feeding
%[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, b_C, r), [1.6  6.6], y0, options);
[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, b_C, r), [1.6 6.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, b_C, r), [1.6 6.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, b_C, r), [1.6 6.6], y0, options);
y = [y1; y2; y];
t = [t1; t2; t];
totalpop = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
infandchron = y(:,3) + y(:,4) + y(:,5) + y(:,6);
carcass = y(:,7);


subplot(3,1,1)
plot(t-0.6, totalpop, 'k-.')
hold on
subplot(3,1,2)
plot(t-0.6, infandchron, 'k-.')
hold on
subplot(3,1,3)
plot(t-0.6, carcass, 'k-.')
hold on


%% Implementation as soon as disease is discovered

y0 = Scen*[0.798, 1.202, 0, 0, 0, 0, 0, 0]; %Initial Conditions

%choose the relevant model formulation with or without feeding
%[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, 0, 0), [0.6 1.6], y0, options);
[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, 0, 0), [0.6 1.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, 0, 0), [0.6 1.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, 0, 0), [0.6 1.6], y0, options);
y0 = y(end,:);
y0(1) = y0(1)-Scen*0.002;
y0(2) = y0(2)-Scen*0.002;
y0(3) = Scen*0.002;
y0(4) = Scen*0.002;
y1 = y;
t1 = t;
%choose the relevant model formulation with or without feeding
%[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, 0, 0), [1.6 1.6+implement], y0, options);
[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, 0, 0), [1.6 1.6+implement], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, 0, 0), [1.6 1.6+implement], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, 0, 0), [1.6 1.6+implement], y0, options);
y0 = y(end,:);
y2 = y;
t2 = t;
%choose the relevant model formulation with or without feeding
%[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 0, b_C, r), [1.6+implement  6.6], y0, options);
[t,y] = ode45(@(t,y) ASFModelEstoniaCNTRL(t, y, 1, b_C, r), [1.6+implement  6.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 0, b_C, r), [1.6+implement  6.6], y0, options);
%[t,y] = ode45(@(t,y) ASFModelSpainCNTRL(t, y, 1, b_C, r), [1.6+implement  6.6], y0, options);
y = [y1; y2; y];
t = [t1; t2; t];
totalpop = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
infandchron = y(:,3) + y(:,4) + y(:,5) + y(:,6);
carcass = y(:,7);


subplot(3,1,1)
plot(t-0.6, totalpop, 'k-')
ylim([0 3])                 %choose y limits according to relevant model simulations
ylabel('Total Density')
%yticks([0 2.5 5])
set(gca,'Xticklabels',[])
set(gca,'Position',[0.12 0.7 0.85 0.25])
str = 'A';
dim = [.9 .9 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,1,2)
plot(t-0.6, infandchron,'k-')
ylim([0 0.015])             %choose y limits according to relevant model simulations
ylabel({'Inf. and';'Chron. Density'})
set(gca,'Xticklabels',[])
set(gca,'Position',[0.12 0.4 0.85 0.25])
str = 'B';
dim = [.9 .6 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');
legend('1 year before', '6 months before','Discovered', 'Location', 'East')

subplot(3,1,3)
plot(t-0.6, carcass, 'k-')
ylim([0 0.03])              %choose y limits according to relevant model simulations
ylabel({'Infected';'Carcass Density'})
%yticks([0 0.01 0.02 0.03])
xlabel('Time')
set(gca,'Position',[0.12 0.1 0.85 0.25])
str = 'C';
dim = [.9 .3 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');
