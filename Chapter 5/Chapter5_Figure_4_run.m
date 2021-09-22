%%% Producing Figure Fig 5.4 of the thesis, Chapter 5

%%% Plotting human density and prevalence levels against a varying carrying
%%% capacity of hosts for the modified Switkes model

%%% Uses function Model_Switkes_wHumans which takes as input the carrying
%%% capacity of hosts, K and the total number of ticks per host, n

close all

Kvec = [2,5:5:100];
hprevlist = zeros(length(Kvec),2);
tprevlist = zeros(length(Kvec),2);
S_S = zeros(length(Kvec),2); 
I_S = zeros(length(Kvec),2); 
R_S = zeros(length(Kvec),2);

for j = 1:2
    for i = 1:length(Kvec)
        K = Kvec(i);
        y0 = [0.006*K, 0.204*K, 0.79*K, 94*K, 6*K, 80, 20, 0];
        [t,y] = ode45(@(t,y) Model_Switkes_wHumans(t,y, K, j*100),[0 200],y0);
        S_S(i,j) = y(end,6);     I_S(i,j) = y(end,7);     R_S(i,j) = y(end,8);
        hprevlist(i,j) = 100*(y(end,7)+y(end,8))/(y(end,6)+y(end,7)+y(end,8));
    end
end

figure
subplot(1,2,1)
set(gca,'Position',[0.08, 0.53, 0.4, 0.44]);
str = 'A';
dim = [.44 .87 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
xlabel('Host Carrying Capacity, K')
yyaxis left
plot(Kvec, S_S(:,1), '--', Kvec, I_S(:,1), '-', Kvec, R_S(:,1), ':', 'LineWidth', 1.5)
ylabel('Human Densities')
ylim([0 100])
yyaxis right
plot(Kvec, hprevlist(:,1), '-','LineWidth', 1.5)
ylim([0 50])
set(gca, 'Yticklabels', [])

subplot(1,2,2)
set(gca,'Position',[0.52, 0.53, 0.4, 0.44]);
str = 'B';
dim = [.88 .87 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
xlabel('Host Carrying Capacity, K')
yyaxis left
plot(Kvec, S_S(:,2), '--', Kvec, I_S(:,2), '-', Kvec, R_S(:,2), ':', 'LineWidth', 1.5)
ylim([0 100])
set(gca, 'Yticklabels',[])
yyaxis right
plot(Kvec, hprevlist(:,2), '-', 'LineWidth', 1.5)
ylim([0 50])
ylabel('Human Sero-prevalence (%)')



