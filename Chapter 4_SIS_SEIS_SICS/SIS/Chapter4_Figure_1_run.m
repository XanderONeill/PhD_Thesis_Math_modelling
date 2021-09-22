%%% Script to produce Figure 4.1 in the Thesis
%%% Requires saved data sets gained from SIS_datasets_runs and also uses
%%% the function DETSIS

close all
K = 1000;
y0 = K*[0.995,0.005];
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'NonNegative',1:2);


figure
load('StochAlpha10')
subplot(2,3,4)
for i = 1:100
    f = find(tvec(i,:) == max(tvec(i,:)));
    plot(tvec(i,1:f), sus(i,1:f), 'k-', tvec(i,1:f), inf(i,1:f), 'b-')
    hold on
    xlabel('Time', 'FontSize', 18)
    ylabel('Population Number','FontSize',18)
    xlim([0 4])
    ylim([0 1100])
    set(gca,'Position', [0.08 0.08 0.27 0.42]);
    str = 'B(i)';
    dim = [.3 .39 .02 .05];
    annotation('textbox', dim, 'String', str,'FontSize',15, 'FitBoxToText', 'off', 'EdgeColor', 'w');
end

load('StochAlpha40')
susoutbreak = zeros(100,length(sus));
infoutbreak = zeros(100,length(inf));
a = zeros(100,1);
for i = 1:100
    f = find(tvec(i,:) == max(tvec(i,:)));
    if min(sus(i,1:f)) < 900
        a(i) = i;
        susoutbreak(i,:) = sus(i,:);
        infoutbreak(i,:) = inf(i,:);
    end
end
a = nonzeros(a);

subplot(2,3,5)
for i = 1:length(a)
    f = find(tvec(a(i),:) == max(tvec(a(i),:)));
    plot(tvec(a(i),1:f), susoutbreak(a(i),1:f), 'k-', tvec(a(i),1:f), infoutbreak(a(i),1:f), 'b-')
    hold on
    xlim([0 40])
    ylim([0 1100])
    xlabel('Time', 'FontSize', 18)
    set(gca,'Yticklabel',[]);
    set(gca,'Position', [0.39 0.08 0.27 0.42]);
    str = 'B(ii)';
    dim = [.61 .39 .02 .05];
    annotation('textbox', dim, 'String', str,'FontSize',15, 'FitBoxToText', 'off','EdgeColor', 'w')%, 'BackgroundColor','w');
end

load('StochAlpha80')
susoutbreak = zeros(100,length(sus));
infoutbreak = zeros(100,length(inf));
a = zeros(100,1);
for i = 1:100
    f = find(tvec(i,:) == max(tvec(i,:)));
    if min(sus(i,1:f)) < 900
        a(i) = i;
        susoutbreak(i,:) = sus(i,:);
        infoutbreak(i,:) = inf(i,:);
    end
end
a = nonzeros(a);

subplot(2,3,6)
for i = 1:length(a)
    f = find(tvec(a(i),:) == max(tvec(a(i),:)));
    plot(tvec(a(i),1:f), susoutbreak(a(i),1:f), 'k-', tvec(a(i),1:f), infoutbreak(a(i),1:f), 'b-')
    hold on
    xlabel('Time', 'FontSize', 18)
    xlim([0 4])
    ylim([0 1100])
    set(gca,'Yticklabel',[]);
    set(gca,'Position', [0.7 0.08 0.27 0.42]);
    str = 'B(iii)';
    dim = [.92 .39 .02 .05];
    annotation('textbox', dim, 'String', str,'FontSize',15, 'FitBoxToText', 'off', 'EdgeColor', 'w');
end

subplot(2,3,1)
[t,y] = ode15s(@(t,y)  DETSIS(t, y, 10), [0 10], y0, options);
sus = y(:,1);inf = y(:,2);
plot(t, sus, 'k-', t, inf, 'b-','LineWidth',1.5)
xlim([0 4])
ylim([0 1100]);
ylabel('Population Density','FontSize',20)
set(gca,'Xticklabel',[]);
set(gca,'Position', [0.08 0.54 0.27 0.42]);
str = 'A(i)';
dim = [.3 .85 .02 .05];
annotation('textbox', dim, 'String', str,'FontSize',15, 'FitBoxToText', 'off', 'EdgeColor', 'w');
    
subplot(2,3,2)
[t,y] = ode15s(@(t,y)  DETSIS(t, y, 40), [0 40], y0, options);
sus = y(:,1);inf = y(:,2);
plot(t, sus, 'k-', t, inf, 'b-','LineWidth',1.5)
xlim([0 40])
ylim([0 1100]);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]);
set(gca,'Position', [0.39 0.54 0.27 0.42]);
str = 'A(ii)';
dim = [.61 .85 .02 .05];
annotation('textbox', dim, 'String', str, 'FontSize',15,'FitBoxToText', 'off', 'EdgeColor', 'w');

subplot(2,3,3)
[t,y] = ode15s(@(t,y)  DETSIS(t, y, 80), [0 10], y0, options);
sus = y(:,1);inf = y(:,2);
plot(t, sus, 'k-', t, inf, 'b-','LineWidth',1.5)
xlim([0 4])
ylim([0 1100]);
set(gca,'Xticklabel',[]);
set(gca,'Yticklabel',[]);
set(gca,'Position', [0.70 0.54 0.27 0.42]);
str = 'A(iii)';
dim = [.92 .85 .02 .05];
annotation('textbox', dim, 'String', str,'FontSize',15, 'FitBoxToText', 'off', 'EdgeColor', 'w');