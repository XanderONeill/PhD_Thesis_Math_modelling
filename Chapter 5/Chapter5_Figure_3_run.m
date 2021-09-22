%%% Producing Figure Fig 5.3 of the thesis, Chapter 5

%%% Plotting tick density and prevalence levels against a varying carrying
%%% capacity of hosts for the modified Switkes model

%%% Uses function Model_Modified_Switkes which takes as input the carrying
%%% capacity of hosts, K and the total number of ticks per host, n

close all

Kvec = [2,5:5:100];
hprevlist = zeros(length(Kvec),2);
tprevlist = zeros(length(Kvec),2);
S_H = zeros(length(Kvec),2); 
I_H = zeros(length(Kvec),2); 
R_H = zeros(length(Kvec),2);

S_T = zeros(length(Kvec),2);
I_T = zeros(length(Kvec),2);

for j = 1:2
    for i = 1:length(Kvec)
        K = Kvec(i);
        y0 = [0.006*K, 0.204*K, 0.79*K, 94*K, 6*K];
        [t,y] = ode45(@(t,y) Model_Modified_Switkes(t,y, K, j*100),[0 200],y0);
        S_H(i,j) = y(end,1);     I_H(i,j) = y(end,2);     R_H(i,j) = y(end,3);
        S_T(i,j) = y(end,4);     I_T(i,j) = y(end,5);
        tprevlist(i,j) = 100*y(end,5)/(y(end,4)+y(end,5));
        hprevlist(i,j) = 100*(y(end,2)+y(end,3))/(y(end,1)+y(end,2)+y(end,3));
    end
end

figure
%%% subplot 1
subplot(2,2,1)
set(gca,'Position',[0.1, 0.57, 0.38, 0.42]);
set(gca, 'Xticklabels',[])
str = 'A(i)';
dim = [.15 .95 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
yyaxis left
plot(Kvec, S_H(:,1), '--', Kvec, I_H(:,1), '-', Kvec, R_H(:,1), ':', 'LineWidth', 1.5)
ylabel('Host Densities')
ylim([0 100])
yyaxis right
plot(Kvec, hprevlist(:,1), 'LineWidth', 1.5)
set(gca, 'Yticklabels',[])
ylim([0 100])

%%% subplot 2
subplot(2,2,2)
set(gca,'Position',[0.53, 0.57, 0.38, 0.42]);
set(gca, 'Xticklabels',[])
str = 'B(i)';
dim = [.61 .95 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
yyaxis left
plot(Kvec, S_H(:,2), '--', Kvec, I_H(:,2), '-', Kvec, R_H(:,2), ':', 'LineWidth', 1.5)
ylim([0 100])
yyaxis right
plot(Kvec, hprevlist(:,2), 'LineWidth', 1.5)
ylabel('Host Prevalence (%)')
ylim([0 100])

%%% subplot 3
subplot(2,2,3)
set(gca,'Position',[0.1, 0.1, 0.38, 0.42]);
str = 'A(ii)';
dim = [.15 .48 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
yyaxis left
plot(Kvec, S_T(:,1), '--', Kvec, I_T(:,1), '-', 'LineWidth', 1.5)
ylabel('Tick Densities')
xlabel('Host Carrying Capacity, K')
yyaxis right
plot(Kvec, tprevlist(:,1), 'LineWidth', 1.5)
ylim([0 100])
set(gca, 'Yticklabels', [])
xlabel('Host Carrying Capacity, K')

%%% subplot 4
subplot(2,2,4)
set(gca,'Position',[0.53, 0.1, 0.38, 0.42]);
str = 'B(ii)';
dim = [.61 .48 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
yyaxis left
plot(Kvec, S_T(:,2), '--', Kvec, I_T(:,2), '-', 'LineWidth', 1.5)
xlabel('Host Carrying Capacity, K')
yyaxis right
plot(Kvec, tprevlist(:,2), 'LineWidth', 1.5)
ylim([0 100])
ylabel('Tick Prevalence (%)')

%plot(2:5:102, susticks, 'b-.', 2:5:102, infticks, 'm-.');
%legend('Susceptible Ticks', 'Infected Ticks','Location','northwest');
%xlabel('Host Carrying Capacity');
%ylabel('Tick Density');
%xlim([0 100])
%ylim([-1 5001])

%plot(2:5:102, hprevlist, 'k-', 2:5:102, tprevlist, 'b-');
%legend('Host Prevalence', 'Tick Prevalence');
%xlabel('Host`s Carrying Capacity')
%ylabel('Prevalence %')
%xlim([0 102])
%ylim([-1 10])