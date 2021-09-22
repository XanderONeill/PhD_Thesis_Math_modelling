%  Chapter 2: Figure 1
%  Model scenarios of Estonia and Spain, under natural conditions and with
%  supplemented feeding

clearvars

options = odeset('Refine',1, 'RelTol',1e-5,'AbsTol',1e-5,'NonNegative',1:8); 
                %Ode solver tolerances and restrictions
            
y0 = [0.796, 1.2, 0.002, 0.002, 0, 0, 0, 0]; %Initial Conditions

%% Estonia model under natural conditions
% running model, plotting total densities, infected densities and ASF
% prevalence

figure              %open figure

[t,y] = ode45(@(t,y) ASFModelEstonia(t, y, 0), [0.6 10.6], y0, options);
totalpop = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
infected = y(:,3) + y(:,4);
chronic = y(:,5) + y(:,6);
prevalence = 100*infected./totalpop;

subplot(3,4,1)                  %total density plot
plot(t-0.6, totalpop, 'k-')
ylim([0 2.5])
ylabel('Total Density')
xticks([0 2 4 6 8 10])
set(gca,'Position', [0.1 0.7 0.18 0.27])
set(gca, 'Xticklabels', []);
str = 'A(i)';
dim = [.21 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');


subplot(3,4,5)                  %infected and chronic density plot
plot(t-0.6, infected, 'k-', t-0.6, chronic, 'k-.');
ylim([0 0.24])
ylabel({'Infected and';'Chronic Density'})
xticks([0 2 4 6 8 10])
set(gca,'Position', [0.1 0.39 0.18 0.27])
set(gca, 'Xticklabels', []);
str = 'A(ii)';
dim = [.21 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,4,9)                 %prevalence plot
plot(t-0.6, prevalence, 'k-');
ylim([0 10]);
yticks([0 2 4 6 8 10]);
ylabel('Prevalence %')
xticks([0 2 4 6 8 10])
xlabel('Time')
set(gca,'Position', [0.1 0.08 0.18 0.27])
str = 'A(iii)';
dim = [.21 .29 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

%% Estonia with supplemented feeding
% running model, plotting total densities, infected densities and ASF
% prevalence

[t,y] = ode45(@(t,y) ASFModelEstonia(t, y, 1), [0.6 10.6], 2*y0, options);
totalpop = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
infected = y(:,3) + y(:,4);
chronic = y(:,5) + y(:,6);
prevalence = 100*infected./totalpop;

subplot(3,4,2)                  %total density plot
plot(t-0.6, totalpop, 'k-')
ylim([0 5])
xticks([0 2 4 6 8 10])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.33 0.7 0.18 0.27])
str = 'B(i)';
dim = [.44 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,4,6)                  %infected and chronic density plot
plot(t-0.6, infected, 'k-', t-0.6, chronic, 'k-.');
ylim([0 0.48])
xticks([0 2 4 6 8 10]);
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.33 0.39 0.18 0.27])
str = 'B(ii)';
dim = [.44 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,4,10)                 %prevalence plot
plot(t-0.6, prevalence, 'k-');
ylim([0 10]);
yticks([0 2 4 6 8 10]);
xlabel('Time')
xticks([0 2 4 6 8 10])
set(gca,'Position', [0.33 0.08 0.18 0.27])
str = 'B(iii)';
dim = [.44 .29 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

%% Spain under natural conditions
% running model, plotting total densities, infected densities and ASF
% prevalence

[t,y] = ode45(@(t,y) ASFModelSpain(t, y, 0), [0.6 10.6], 2.5*y0, options);
totalpop = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
infected = y(:,3) + y(:,4);
chronic = y(:,5) + y(:,6);
prevalence = 100*infected./totalpop;

subplot(3,4,3)
plot(t-0.6, totalpop, 'k-')
ylim([0 6.25]);
yticks([0 1 2 3 4 5 6]);
xticks([0 2 4 6 8 10]);
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.56 0.7 0.18 0.27]);
str = 'C(i)';
dim = [.67 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,4,7)
plot(t-0.6, infected, 'k-', t-0.6, chronic, 'k-.');
ylim([0 0.6])
xticks([0 2 4 6 8 10]);
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.56 0.39 0.18 0.27]);
str = 'C(ii)';
dim = [.67 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,4,11)
plot(t-0.6, prevalence, 'k-');
ylim([0 10]);
yticks([0 2 4 6 8 10]);
xlabel('Time');
xticks([0 2 4 6 8 10]);
set(gca,'Position', [0.56 0.08 0.18 0.27]);
str = 'C(iii)';
dim = [.67 .29 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

%% Spain with supplemented feeding
% running model, plotting total densities, infected densities and ASF
% prevalence

[t,y] = ode45(@(t,y) ASFModelSpain(t, y, 1), [0.6 10.6], 5*y0, options);
totalpop = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
infected = y(:,3) + y(:,4);
chronic = y(:,5) + y(:,6);
prevalence = 100*infected./totalpop;

subplot(3,4,4)
plot(t-0.6, totalpop, 'k-')
ylim([0 12.5]);
yticks([0 2 4 6 8 10 12]);
xticks([0 2 4 6 8 10]);
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.79 0.70 0.18 0.27]);
str = 'D(i)';
dim = [.90 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,4,8)
plot(t-0.6, infected, 'k-', t-0.6, chronic, 'k-.');
ylim([0 1.2])
xticks([0 2 4 6 8 10]);
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.79 0.39 0.18 0.27])
str = 'D(ii)';
dim = [.90 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,4,12)
plot(t-0.6, prevalence, 'k-');
ylim([0 10]);
yticks([0 2 4 6 8 10]);
xlabel('Time');
xticks([0 2 4 6 8 10]);
set(gca,'Position', [0.79 0.08 0.18 0.27]);
str = 'D(iii)';
dim = [.90 .29 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');
