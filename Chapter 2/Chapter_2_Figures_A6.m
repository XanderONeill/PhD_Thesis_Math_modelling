%  Chapter 2: Figure A.6
%  Model results for Estonia under natural conditions with different values
%  of survivor reversion rates (along with parameter values that satisfy
%  the epidemiological criteria)

clearvars
options = odeset('Refine',1, 'RelTol',1e-8,'AbsTol',1e-8,'NonNegative',1:8); 
                %Ode solver tolerances and restrictions
            
y0 = [0.796, 1.2, 0.002, 0.002, 0, 0, 0, 0]; %Initial Conditions


%% Varying survivor reversion rate kappa

figure              %open figure
%%% running the model and determining densities

% parameter vectors that satisfy the epidemiological criteria for different
% values of kappa
kappa_vec =     [12/2, 12/3, 12/4, 12/6, 12/9];
beta_F_vec =    [40, 50, 58, 63, 67.5];
beta_E_vec =    [4, 3.1, 2.1, 2, 1.7];
rho_vec =       [0.61, 0.7, 0.8, 0.85, 0.89];
for i = 1:5
    kappa = kappa_vec(i);
    beta_F = beta_F_vec(i);
    beta_E = beta_E_vec(i);
    rho = rho_vec(i);
    [t,y] = ode45(@(t,y) ASFModelEstoniaFigA6(t, y, 0, beta_F, beta_E, rho, kappa), [0.6 5.6], y0, options);
    totalpop = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
    infected = y(:,3) + y(:,4);
    chronic = y(:,5) + y(:,6);
    prevalence = 100*infected./totalpop;

    subplot(3,5,i)                  %total density plot
    plot(t-0.6, totalpop, 'k-')

    subplot(3,5,5+i)                  %infected and chronic density plot
    plot(t-0.6, infected, 'k-', t-0.6, chronic, 'k-.');
    
    subplot(3,5,10+i)                 %prevalence plot
    plot(t-0.6, prevalence, 'k-');

end

%% Plots: [1]
subplot(3,5,1)
ylim([0 2.5])
ylabel('Total Density')
xlim([0 5])
xticks([0 2 4])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.1 0.7 0.14 0.25])
str = 'A(i)';
dim = [.17 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,5,6)
ylim([0 0.25])
yticks([0 0.05 0.1 0.15 0.2])
ylabel({'Inf. and';'Chron. Density'})
xlim([0 5])
xticks([0 2 4])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.1 0.4 0.14 0.25])
str = 'A(ii)';
dim = [.17 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,5,11)
ylim([0 6.5])
yticks([0 2 4 6])
ylabel({'Prevalence (%)'})
xlim([0 5])
xlabel('Time (years)')
xticks([0 2 4])
set(gca,'Position', [0.1 0.1 0.14 0.25])
str = 'A(iii)';
dim = [.17 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

%% Plots: [2]
subplot(3,5,2)
ylim([0 2.5])
set(gca, 'Yticklabels', []);
xlim([0 5])
xticks([0 2 4])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.27 0.7 0.14 0.25])
str = 'B(i)';
dim = [.34 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,5,7)
ylim([0 0.25])
yticks([0 0.05 0.1 0.15 0.2])
set(gca, 'Yticklabels', []);
xlim([0 5])
xticks([0 2 4])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.27 0.4 0.14 0.25])
str = 'B(ii)';
dim = [.34 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');
    
subplot(3,5,12)
ylim([0 6.5])
yticks([0 2 4 6])
set(gca, 'Yticklabels', []);
xlim([0 5])
xlabel('Time (years)')
xticks([0 2 4])
set(gca,'Position', [0.27 0.1 0.14 0.25])
str = 'B(iii)';
dim = [.34 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

%% Plots [3]
subplot(3,5,3)
ylim([0 2.5])
set(gca, 'Yticklabels', []);
xlim([0 5])
xticks([0 2 4])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.44 0.7 0.14 0.25])
str = 'C(i)';
dim = [.51 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');   

subplot(3,5,8)
ylim([0 0.25])
yticks([0 0.05 0.1 0.15 0.2])
set(gca, 'Yticklabels', []);
xlim([0 5])
xticks([0 2 4])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.44 0.4 0.14 0.25])
str = 'C(ii)';
dim = [.51 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,5,13)
ylim([0 6.5])
yticks([0 2 4 6])
set(gca, 'Yticklabels', []);
xlim([0 5])
xlabel('Time (years)')
xticks([0 2 4])
set(gca,'Position', [0.44 0.1 0.14 0.25])
str = 'C(iii)';
dim = [.51 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');
   
%% Plots [4]
subplot(3,5,4)
ylim([0 2.5])
set(gca, 'Yticklabels', []);
xlim([0 5])
xticks([0 2 4])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.61 0.7 0.14 0.25])
str = 'D(i)';
dim = [.68 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none'); 

subplot(3,5,9)
ylim([0 0.25])
yticks([0 0.05 0.1 0.15 0.2])
set(gca, 'Yticklabels', []);
xlim([0 5])
xticks([0 2 4])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.61 0.4 0.14 0.25])
str = 'D(ii)';
dim = [.68 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,5,14)
ylim([0 6.5])
yticks([0 2 4 6])
set(gca, 'Yticklabels', []);
xlim([0 5])
xlabel('Time (years)')
xticks([0 2 4])
set(gca,'Position', [0.61 0.1 0.14 0.25])
str = 'D(iii)';
dim = [.68 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

%% Plots [5]
subplot(3,5,5)
ylim([0 2.5])
set(gca, 'Yticklabels', []);
xlim([0 5])
xticks([0 2 4])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.78 0.7 0.14 0.25])
str = 'E(i)';
dim = [.85 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none'); 

subplot(3,5,10)
ylim([0 0.25])
yticks([0 0.05 0.1 0.15 0.2])
set(gca, 'Yticklabels', []);
xlim([0 5])
xticks([0 2 4])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.78 0.4 0.14 0.25])
str = 'E(ii)';
dim = [.85 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,5,15)
ylim([0 6.5])
yticks([0 2 4 6])
set(gca, 'Yticklabels', []);
xlim([0 5])
xlabel('Time (years)')
xticks([0 2 4])
set(gca,'Position', [0.78 0.1 0.14 0.25])
str = 'E(iii)';
dim = [.85 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');