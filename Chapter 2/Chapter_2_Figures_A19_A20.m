%%% Exploring "Effect" of Vultures by increasing and decreasing degradation
%%% rate of infected carcasses.

%%% Uses model ASFModelSpainDEG, and can be used to produce figures:
%%% Fig A.19 and Fig A.20 in Chapter 2 of the thesis.

%%% There are two sections to this script, the first is to produce
%%% Fig A.19 and the second for Fig A.20

clearvars

options = odeset('Refine',1, 'RelTol',1e-8,'AbsTol',1e-9, 'NonNegative',1:8);
        %Ode solver tolerances and restrictions

y0 = 2.5*[0.796, 1.2, 0.002, 0.002, 0, 0, 0, 0]; 
        %Initial conditions for the system representative of Spain under
        %natural conditions. Adapt accordingly for other scenarios.

%% Fig A.19
%%% 1 day -> 1/2 week -> 1 week (default) degradation

%%% 1 day degradation
[t,y] = ode45(@(t,y) ASFModelSpainDEG(t, y, 0, 1/7), [0.6 5.6], y0, options);
t4 = t;
totalpop4 = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
prevalence4 = 100*(y(:,3) + y(:,4))./totalpop4;
infected4 = y(:,3) + y(:,4);
chronic4 = y(:,5) + y(:,6);

figure
subplot(4,3,1)
plot(t4-0.6,totalpop4,'k-')
ylim([0 9])
ylabel('Tot. Pop.')
yticks([2 4 6 8])
xlim([0 5])
set(gca,'XTickLabels',[]);
set(gca,'Position',[0.1 0.77 0.27 0.2])
str = 'A(i)';
dim = [.3 .94 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,4)
plot(t4-0.6,infected4,'k-', t4-0.6,chronic4,'k-.')
ylim([0 0.15])
ylabel({'Inf. and';'Chron. Pop'})
yticks([0.05 0.1 0.15])
xlim([0 5])
set(gca,'XTickLabels',[]);
set(gca,'Position',[0.1 0.54 0.27 0.2])
str = 'A(ii)';
dim = [.3 .71 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,7)
plot(t4-0.6,y(:,7),'k-')
ylim([0 0.06])
ylabel({'Inf. Carcass';'Density'})
yticks([0.02 0.04 0.06])
xlim([0 5])
set(gca,'XTickLabels',[]);
set(gca,'Position',[0.1 0.31 0.27 0.2])
str = 'A(iii)';
dim = [.3 .48 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,10)
plot(t4-0.6,prevalence4,'k-')
ylim([0 1.5])
yticks([0.5 1 1.5])
ylabel('Prev. (%)')
xlabel('Time')
xlim([0 5])
set(gca,'Position',[0.1 0.08 0.27 0.2])
str = 'A(iv)';
dim = [.3 .25 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');


%%% 0.5 week degradation
[t,y] = ode45(@(t,y) ASFModelSpainDEG(t, y, 0, 0.5), [0.6 5.6], y0, options);
t5 = t;
totalpop5 = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
prevalence5 = 100*(y(:,3) + y(:,4))./totalpop5;
infected5 = y(:,3) + y(:,4);
chronic5 = y(:,5) + y(:,6);

subplot(4,3,2)
plot(t5-0.6,totalpop5,'k-')
ylim([0 9])
yticks([2 4 6 8])
set(gca,'YTickLabels',[]);
xlim([0 5])
set(gca,'XTickLabels',[]);
set(gca,'Position',[0.4 0.77 0.27 0.2])
str = 'B(i)';
dim = [.6 .94 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,5)
plot(t5-0.6,infected5,'k-', t5-0.6,chronic5,'k-.')
ylim([0 0.15])
yticks([0.05 0.1 0.15])
set(gca,'YTickLabels',[]);
xlim([0 5])
set(gca,'XTickLabels',[]);
legend('Inf.','Chron.','Location','NorthWest')
set(gca,'Position',[0.4 0.54 0.27 0.2])
str = 'B(ii)';
dim = [.6 .71 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,8)
plot(t5-0.6,y(:,7),'k-')
ylim([0 0.06])
yticks([0.02 0.04 0.06])
set(gca,'YTickLabels',[]);
xlim([0 5])
set(gca,'XTickLabels',[]);
set(gca,'Position',[0.4 0.31 0.27 0.2])
str = 'B(iii)';
dim = [.6 .48 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,11)
plot(t5-0.6,prevalence5,'k-')
ylim([0 1.5])
yticks([0.5 1 1.5])
set(gca,'YTickLabels',[]);
xlim([0 5])
xlabel('Time')
set(gca,'Position',[0.4 0.08 0.27 0.2])
str = 'B(iv)';
dim = [.6 .25 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');

%%% 1 week degradation
[t,y] = ode45(@(t,y) ASFModelSpainDEG(t, y, 0, 1), [0.6 5.6], y0, options);
t6 = t;
totalpop6 = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
prevalence6 = 100*(y(:,3) + y(:,4))./totalpop6;
infected6 = y(:,3) + y(:,4);
chronic6 = y(:,5) + y(:,6);

subplot(4,3,3)
plot(t6-0.6,totalpop6,'k-')
ylim([0 9])
yticks([2 4 6 8])
set(gca,'YTickLabels',[]);
xlim([0 5])
set(gca,'XTickLabels',[]);
set(gca,'Position',[0.7 0.77 0.27 0.2])
str = 'C(i)';
dim = [.9 .94 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,6)
plot(t6-0.6,infected6,'k-', t6-0.6,chronic6,'k-.')
ylim([0 0.15])
yticks([0.05 0.1 0.15])
set(gca,'YTickLabels',[]);
xlim([0 5])
set(gca,'XTickLabels',[]);
set(gca,'Position',[0.7 0.54 0.27 0.2])
str = 'C(ii)';
dim = [.9 .71 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,9)
plot(t6-0.6,y(:,7),'k-')
ylim([0 0.06])
yticks([0.02 0.04 0.06])
set(gca,'YTickLabels',[]);
xlim([0 5])
set(gca,'XTickLabels',[]);
set(gca,'Position',[0.7 0.31 0.27 0.2])
str = 'C(iii)';
dim = [.9 .48 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,12)
plot(t6-0.6,prevalence6,'k-')
ylim([0 1.5])
yticks([0.5 1 1.5])
set(gca,'YTickLabels',[]);
xlim([0 5])
xlabel('Time')
set(gca,'Position',[0.7 0.08 0.27 0.2])
str = 'C(iv)';
dim = [.9 .25 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
        
%% Fig A.20
% 1 week (default) -> 2 weeks -> 4 weeks degradation

%%% 1 week degradation
[t,y] = ode45(@(t,y) ASFModelSpainDEG(t, y, 0, 1), [0.6 5.6], y0, options);
t1 = t;
totalpop1 = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
prevalence1 = 100*(y(:,3) + y(:,4))./totalpop1;
infected1 = y(:,3) + y(:,4);
chronic1 = y(:,5) + y(:,6);

figure
subplot(4,3,1)
plot(t1-0.6,totalpop1,'k-')
ylim([0 6])
ylabel('Tot. Pop.')
xlim([0 5])
set(gca,'XTickLabels',[]);
set(gca,'Position',[0.1 0.77 0.27 0.2])
str = 'A(i)';
dim = [.3 .94 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,4)
plot(t1-0.6,infected1,'k-', t1-0.6,chronic1,'k-.')
xlim([0 5])
ylim([0 0.3])
ylabel({'Inf. and';'Chron. Pop'})
set(gca,'Position',[0.1 0.54 0.27 0.2])
set(gca,'XTickLabels',[]);
str = 'A(ii)';
dim = [.3 .71 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,7)
plot(t1-0.6,y(:,7),'k-')
xlim([0 5])
xlim([0 5])
ylim([0 0.60])
ylabel({'Inf. Carcass';'Density'})
set(gca,'XTickLabels',[]);
set(gca,'Position',[0.1 0.31 0.27 0.2])
str = 'A(iii)';
dim = [.3 .48 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,10)
plot(t1-0.6,prevalence1,'k-')
xlim([0 5])
xlim([0 5])
ylim([0 6])
ylabel('Prev. (%)')
xlabel('Time')
set(gca,'Position',[0.1 0.08 0.27 0.2])
str = 'A(iv)';
dim = [.3 .25 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');

%%% 2 week degradation
[t,y] = ode45(@(t,y) ASFModelSpainDEG(t, y, 0, 2), [0.6 5.6], y0, options);
t2 = t;
totalpop2 = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
prevalence2 = 100*(y(:,3) + y(:,4))./totalpop2;
infected2 = y(:,3) + y(:,4);
chronic2 = y(:,5) + y(:,6);

subplot(4,3,2)
plot(t2-0.6,totalpop2,'k-')
xlim([0 5])
ylim([0 6])
set(gca,'Position',[0.4 0.77 0.27 0.2])
set(gca,'XTickLabels',[]);
set(gca,'YTickLabels',[]);
str = 'B(i)';
dim = [.6 .94 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,5)
plot(t2-0.6,infected2,'k-', t2-0.6,chronic2,'k-.')
xlim([0 5])
ylim([0 0.3])
legend('Inf.','Chron.', 'Location', 'East')
set(gca,'Position',[0.4 0.54 0.27 0.2])
set(gca,'XTickLabels',[]);
set(gca,'YTickLabels',[]);
str = 'B(ii)';
dim = [.6 .71 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,8)
plot(t2-0.6,y(:,7),'k-')
xlim([0 5])
ylim([0 0.6])
set(gca,'XTickLabels',[]);
set(gca,'YTickLabels',[]);
set(gca,'Position',[0.4 0.31 0.27 0.2])
str = 'B(iii)';
dim = [.6 .48 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,11)
plot(t2-0.6,prevalence2,'k-')
xlim([0 5])
ylim([0 6])
xlabel('Time')
set(gca,'YTickLabels',[]);
set(gca,'Position',[0.4 0.08 0.27 0.2])
str = 'B(iv)';
dim = [.6 .25 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');

%%% 4 week degradation
[t,y] = ode45(@(t,y) ASFModelSpainDEG(t, y, 0, 4), [0.6 5.6], y0, options);
t3 = t;
totalpop3 = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
prevalence3 = 100*(y(:,3) + y(:,4))./totalpop3;
infected3 = y(:,3) + y(:,4);
chronic3 = y(:,5) + y(:,6);

%%%
subplot(4,3,3)
plot(t3-0.6,totalpop3,'k-')
xlim([0 5])
ylim([0 6])
set(gca,'Position',[0.7 0.77 0.27 0.2])
set(gca,'XTickLabels',[]);
set(gca,'YTickLabels',[]);
str = 'C(i)';
dim = [.9 .94 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,6)
plot(t3-0.6,infected3,'k-', t-0.6,chronic3,'k-.')
xlim([0 5])
ylim([0 0.3])
set(gca,'Position',[0.7 0.54 0.27 0.2])
set(gca,'XTickLabels',[]);
set(gca,'YTickLabels',[]);
str = 'C(ii)';
dim = [.9 .71 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,9)
plot(t3-0.6,y(:,7),'k-')
xlim([0 5])
ylim([0 0.60])
set(gca,'XTickLabels',[]);
set(gca,'YTickLabels',[]);
set(gca,'Position',[0.7 0.31 0.27 0.2])
str = 'C(iii)';
dim = [.9 .48 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
%%%
subplot(4,3,12)
plot(t3-0.6,prevalence3,'k-')
xlim([0 5])
ylim([0 6])
xlabel('Time')
set(gca,'YTickLabels',[]);
set(gca,'Position',[0.7 0.08 0.27 0.2])
str = 'C(iv)';
dim = [.9 .25 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
