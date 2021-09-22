%%% Producing Figure Fig 5.2 of the thesis, Chapter 5

%%% Plotting prevalence levels against a varying number of ticks that can
%%% be supported by hosts

%%% Uses function Model_Modified_Switkes which takes as input the carrying
%%% capacity of hosts, K and the total number of ticks per host, n

nvec = 0:200;
hprevlist = zeros(1,length(nvec));
tprevlist = zeros(1,length(nvec));

for j = 1:length(nvec)
   y0 = [0, 2.04, 7.88, 0.938*nvec(j), 0.062*nvec(j)];
   [t,y] = ode45(@(t,y) Model_Modified_Switkes(t,y, 10, nvec(j)),[0 100],y0);
   tprevlist(1,j) = (100*y(end,5)./(y(end,4)+y(end,5)));
   hprevlist(1,j) = (100*(y(end,2)+y(end,3))./(y(end,1)+y(end,2)+y(end,3)));
end


figure
yyaxis left
plot(0:200, hprevlist(1,:), '-', 'LineWidth', 1.5);
set(gca,'Position',[0.1, 0.1, 0.7, 0.76]);
xlabel('Average number of ticks that can be supported by a host, \it{n}')
ylabel('Prevalence %')
xlim([0 200])
ylim([0 100])
yyaxis right
plot(0:200, tprevlist(1,:), '-', 'LineWidth', 1.5);
set(gca, 'Yticklabels',[]);
xlim([0 200])
xticks([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
ylim([0 100])

