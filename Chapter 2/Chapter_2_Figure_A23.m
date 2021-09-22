%%% Script used to produce Figure Fig A.23 in Chapter 2 of the thesis

%%% Here, we show how the epidemiological criteria for an ASF outbreak
%%% cannot be matched for a varying carcass degradation rate without a
%%% survivor class

%%% Running Simulations representative of Estonia under natural conditions
%%% with transmission coefficients to show the possible outcomes seen for
%%% different degradation rates

clearvars
options = odeset('Refine',1, 'RelTol',1e-8,'AbsTol',1e-9, 'NonNegative',1:7);
            % ode model tolerances and restrictions
y0 = [0.796, 1.2, 0.002, 0.002, 0, 0, 0];
            %initial conditions for Estonia under natural conditions
%% General Code

[t,y] = ode45(@(t,y) ASFModelEstoniaNOSURV(t, y, 0, 40, 2, 25), [0.6 5.6], y0, options);
totalpop = y(:,1) + y(:,2)+ y(:,3) + y(:,4) + y(:,5) + y(:,6);
infected = y(:,3) + y(:,4);
prevalence = 100*infected./totalpop;

figure
subplot(3,2,1)
plot(t-0.6, totalpop, 'k-')
ylabel('Total Density')
ylim([0 3])
set(gca,'XTickLabels',[])
set(gca,'Position',[0.1 0.68 0.42 0.28])
str = 'A(i)';
dim = [.45 .89 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');

subplot(3,2,3)
plot(t-0.6, infected, 'k-')
ylabel('Infected Density')
ylim([0 0.04])
set(gca,'XTickLabels',[])
set(gca, 'Position', [0.1 0.38 0.42 0.28])
str = 'A(ii)';
dim = [.45 .59 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');

subplot(3,2,5)
plot(t-0.6, prevalence, 'k-')
ylabel('Prevalence (%)')
ylim([0 4.5])
xlabel('Time')
set(gca, 'Position', [0.1 0.08 0.42 0.28])
str = 'A(iii)';
dim = [.45 .29 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');


clearvars   %resetting variables to run model for new parameters
options = odeset('Refine',1, 'RelTol',1e-8,'AbsTol',1e-9, 'NonNegative',1:7);

y0 = [0.796, 1.2, 0.002, 0.002, 0, 0, 0];
[t,y] = ode45(@(t,y) ASFModelEstoniaNOSURV(t, y, 0, 42, 1, 40), [0.6 5.6], y0, options);        
totalpop = y(:,1) + y(:,2)+ y(:,3) + y(:,4) + y(:,5) + y(:,6);
infected = y(:,3) + y(:,4);
chronic = y(:,5) + y(:,6);
prevalence = 100*infected./totalpop;

subplot(3,2,2)
plot(t-0.6, totalpop, 'k-')
ylim([0 3])
set(gca,'YTickLabels',[])
set(gca,'XTickLabels',[])
set(gca,'Position',[0.55 0.68 0.42 0.28])
str = 'B(i)';
dim = [.9 .89 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');

subplot(3,2,4)
plot(t-0.6, infected,'k-')
ylim([0 0.04])
set(gca,'YTickLabels',[])
set(gca,'XTickLabels',[])
set(gca, 'Position', [0.55 0.38 0.42 0.28])
str = 'B(ii)';
dim = [.9 .59 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');

subplot(3,2,6)
plot(t-0.6, prevalence, 'k-')
ylim([0 4.5])
set(gca,'YTickLabels',[])
xlabel('Time')
set(gca, 'Position', [0.55 0.08 0.42 0.28])
str = 'B(iii)';
dim = [.9 .29 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');