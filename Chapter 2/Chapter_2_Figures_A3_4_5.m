%  Chapter 2: Figure A.3 and A.4 and A.5
%  Model results for Estonia under natural conditions with different values
%  of transmission coefficients and rho

clearvars
options = odeset('Refine',1, 'RelTol',1e-8,'AbsTol',1e-8,'NonNegative',1:8); 
                %Ode solver tolerances and restrictions
            
y0 = [0.796, 1.2, 0.002, 0.002, 0, 0, 0, 0]; %Initial Conditions

%% Varying environmental transmission coefficient beta_E

figure              %open figure
%%% running the model and determining densities
beta_E_vec = [0, 2, 6];
for i = 1:3
    beta_E = beta_E_vec(i);
    [t,y] = ode45(@(t,y) ASFModelEstoniaCRTR(t, y, 0, 63, beta_E, 0.85), [0.6 10.6], y0, options);
    totalpop = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
    infected = y(:,3) + y(:,4);
    chronic = y(:,5) + y(:,6);
    prevalence = 100*infected./totalpop;

    subplot(3,3,i)                  %total density plot
    plot(t-0.6, totalpop, 'k-')

    subplot(3,3,3+i)                  %infected and chronic density plot
    plot(t-0.6, infected, 'k-', t-0.6, chronic, 'k-.');
    
    subplot(3,3,6+i)                 %prevalence plot
    plot(t-0.6, prevalence, 'k-');

end
subplot(3,3,1)
ylim([0 3])
ylabel('Total Density')
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.1 0.7 0.25 0.25])
str = 'A(i)';
dim = [.29 .73 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,2)
ylim([0 3])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.4 0.7 0.25 0.25])
str = 'B(i)';
dim = [.59 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');
    
subplot(3,3,3)
ylim([0 3])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.7 0.7 0.25 0.25])
str = 'C(i)';
dim = [.89 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');   
       
subplot(3,3,4)
ylim([0 0.2])
ylabel({'Inf. and Chron. Dens.'})
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.1 0.4 0.25 0.25])
str = 'A(ii)';
dim = [.29 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,5)
ylim([0 0.2])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.4 0.4 0.25 0.25])
str = 'B(ii)';
dim = [.59 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,6)
ylim([0 0.2])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.7 0.4 0.25 0.25])
str = 'C(ii)';
dim = [.89 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,7)
ylim([0 11])
yticks([0 2 4 6 8 10])
ylabel({'Prevalence (%)'})
xlim([0 8])
xlabel('Time (years)')
xticks([0 2 4 6 8])
set(gca,'Position', [0.1 0.1 0.25 0.25])
str = 'A(iii)';
dim = [.28 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,8)
ylim([0 11])
yticks([0 2 4 6 8 10])
set(gca, 'Yticklabels', []);
xlim([0 8])
xlabel('Time (years)')
xticks([0 2 4 6 8])
set(gca,'Position', [0.4 0.1 0.25 0.25])
str = 'B(iii)';
dim = [.58 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,9)
ylim([0 11])
yticks([0 2 4 6 8 10])
set(gca, 'Yticklabels', []);
xlim([0 8])
xlabel('Time (years)')
xticks([0 2 4 6 8])
set(gca,'Position', [0.7 0.1 0.25 0.25])
str = 'C(iii)';
dim = [.88 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

%% Varying frequency-dependent transmission coefficient beta_F

figure              %open figure
%%% running the model and determining densities
beta_F_vec = [50, 63, 76];
for i = 1:3
    beta_F = beta_F_vec(i);
    [t,y] = ode45(@(t,y) ASFModelEstoniaCRTR(t, y, 0, beta_F, 2, 0.85), [0.6 10.6], y0, options);
    totalpop = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
    infected = y(:,3) + y(:,4);
    chronic = y(:,5) + y(:,6);
    prevalence = 100*infected./totalpop;

    subplot(3,3,i)                  %total density plot
    plot(t-0.6, totalpop, 'k-')

    subplot(3,3,3+i)                  %infected and chronic density plot
    plot(t-0.6, infected, 'k-', t-0.6, chronic, 'k-.');
    
    subplot(3,3,6+i)                 %prevalence plot
    plot(t-0.6, prevalence, 'k-');

end
subplot(3,3,1)
ylim([0 3])
ylabel('Total Density')
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.1 0.7 0.25 0.25])
str = 'A(i)';
dim = [.29 .73 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,2)
ylim([0 3])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.4 0.7 0.25 0.25])
str = 'B(i)';
dim = [.59 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');
    
subplot(3,3,3)
ylim([0 3])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.7 0.7 0.25 0.25])
str = 'C(i)';
dim = [.89 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');   
       
subplot(3,3,4)
ylim([0 0.2])
ylabel({'Inf. and Chron. Dens.'})
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.1 0.4 0.25 0.25])
str = 'A(ii)';
dim = [.29 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,5)
ylim([0 0.2])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.4 0.4 0.25 0.25])
str = 'B(ii)';
dim = [.59 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,6)
ylim([0 0.2])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.7 0.4 0.25 0.25])
str = 'C(ii)';
dim = [.89 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,7)
ylim([0 11])
yticks([0 2 4 6 8 10])
ylabel({'Prevalence (%)'})
xlim([0 8])
xlabel('Time (years)')
xticks([0 2 4 6 8])
set(gca,'Position', [0.1 0.1 0.25 0.25])
str = 'A(iii)';
dim = [.28 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,8)
ylim([0 11])
yticks([0 2 4 6 8 10])
set(gca, 'Yticklabels', []);
xlim([0 8])
xlabel('Time (years)')
xticks([0 2 4 6 8])
set(gca,'Position', [0.4 0.1 0.25 0.25])
str = 'B(iii)';
dim = [.58 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,9)
ylim([0 11])
yticks([0 2 4 6 8 10])
set(gca, 'Yticklabels', []);
xlim([0 8])
xlabel('Time (years)')
xticks([0 2 4 6 8])
set(gca,'Position', [0.7 0.1 0.25 0.25])
str = 'C(iii)';
dim = [.88 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');


%% Varying rho

figure              %open figure
%%% running the model and determining densities
rho_vec = [0.7, 0.85, 1];
for i = 1:3
    rho = rho_vec(i);
    [t,y] = ode45(@(t,y) ASFModelEstoniaCRTR(t, y, 0, 63, 2, rho), [0.6 10.6], y0, options);
    totalpop = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
    infected = y(:,3) + y(:,4);
    chronic = y(:,5) + y(:,6);
    prevalence = 100*infected./totalpop;

    subplot(3,3,i)                  %total density plot
    plot(t-0.6, totalpop, 'k-')

    subplot(3,3,3+i)                  %infected and chronic density plot
    plot(t-0.6, infected, 'k-', t-0.6, chronic, 'k-.');
    
    subplot(3,3,6+i)                 %prevalence plot
    plot(t-0.6, prevalence, 'k-');

end
subplot(3,3,1)
ylim([0 3])
ylabel('Total Density')
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.1 0.7 0.25 0.25])
str = 'A(i)';
dim = [.29 .73 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,2)
ylim([0 3])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.4 0.7 0.25 0.25])
str = 'B(i)';
dim = [.59 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');
    
subplot(3,3,3)
ylim([0 3])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.7 0.7 0.25 0.25])
str = 'C(i)';
dim = [.89 .91 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');   
       
subplot(3,3,4)
ylim([0 0.2])
ylabel({'Inf. and Chron. Dens.'})
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.1 0.4 0.25 0.25])
str = 'A(ii)';
dim = [.29 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,5)
ylim([0 0.2])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.4 0.4 0.25 0.25])
str = 'B(ii)';
dim = [.59 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,6)
ylim([0 0.2])
set(gca, 'Yticklabels', []);
xlim([0 8])
xticks([0 2 4 6 8])
set(gca, 'Xticklabels', []);
set(gca,'Position', [0.7 0.4 0.25 0.25])
str = 'C(ii)';
dim = [.89 .61 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,7)
ylim([0 11])
yticks([0 2 4 6 8 10])
ylabel({'Prevalence (%)'})
xlim([0 8])
xlabel('Time (years)')
xticks([0 2 4 6 8])
set(gca,'Position', [0.1 0.1 0.25 0.25])
str = 'A(iii)';
dim = [.28 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,8)
ylim([0 11])
yticks([0 2 4 6 8 10])
set(gca, 'Yticklabels', []);
xlim([0 8])
xlabel('Time (years)')
xticks([0 2 4 6 8])
set(gca,'Position', [0.4 0.1 0.25 0.25])
str = 'B(iii)';
dim = [.58 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');

subplot(3,3,9)
ylim([0 11])
yticks([0 2 4 6 8 10])
set(gca, 'Yticklabels', []);
xlim([0 8])
xlabel('Time (years)')
xticks([0 2 4 6 8])
set(gca,'Position', [0.7 0.1 0.25 0.25])
str = 'C(iii)';
dim = [.88 .31 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on','EdgeColor', 'none');