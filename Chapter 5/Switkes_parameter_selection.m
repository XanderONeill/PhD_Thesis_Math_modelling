%%% Parameter sensitivity for the modified Switkes model to get the correct
%%% transmission parameters that match with the data (see citations in
%%% Chapter 5).

clearvars
close all

%% Finding the Correct Beta and Gamma Values

y0 = [4, 2, 4, 600, 400];
betaHVvec = 0:0.0001:0.02;
gammavec = 0:0.1:1;
tprevlist=zeros(length(betaHVvec),length(gammavec));
hprevlist=zeros(length(betaHVvec),length(gammavec));
for i = 1:length(betaHVvec)
    for j = 1:length(gammavec)
        betaHV = betaHVvec(i);
        gamma = gammavec(j);
        [t,y] = ode45(@(t,y) Model_Modified_Switkes(t,y,10,betaHV, gamma),[0 100],y0);
        tickprevalence = (100*y(end,5)./(y(end,4)+y(end,5)));
        hostprevalence = (100*(y(end,2)+y(end,3))./(y(end,1)+y(end,2)+y(end,3)));
        tprevlist(i,j) = tickprevalence;
        hprevlist(i,j) = hostprevalence;
    end
end

find(tprevlist<5 & tprevlist>4 & hprevlist<20.5 & hprevlist>19.5)

%list = [0.0101:0.0001:0.013];
%plot(list, tprevlist, 'k-', list, hprevlist, 'b-');
%leg = legend('tick prevalence', 'host prevalence');
%set(leg,'visible','on')
%xlabel('betaHV coefficient');
%ylabel('Prevalence')
%hold on
%h=18; h1=21; h2=2; h3=7;
%plot([0.01 0.013], [h h], 'r-', [0.01 .013], [h1 h1], 'r-', [0.01 .013], [h2 h2], 'y-', [0.01 .013], [h3 h3], 'y-');
%hold off


