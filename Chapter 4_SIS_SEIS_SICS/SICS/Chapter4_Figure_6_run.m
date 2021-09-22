%%% Producing Figure 4.6 in the thesis, chapter 4
%%% Pathogen extinction time plots with error bars, for varying model
%%% parameters

%%% This script is used in conjuction with the similarly named script found
%%% for SIS.

%%% Run the SIS script followed by this

z = 1.96;    %95% confidence interval
bF = 0:10:60;
bD = (0:6)/20;
gam = 5:5:80;
kap = 10:10:160;
eps = 0:0.05:1;
T = 500;



%%% 1-4 beta_F, 5-8 beta_D, 9-12 gamma, 13-16 kappa, 17-18 epsilon
%%% ODDS: epsilon = 0.5, c = 0
%%% EVEN: epsilon = 0, c = 0.5
%% #1                           Beta_F
load('SICSc0Eps05betaFAlpha80it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2)); 
cimin = zeros(1, sz(2));
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];
   
    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,6)
if isempty(g) == 0
    errorbar(bF(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 60])
xticks([0 20 40 60])
xlabel('\beta_F','FontSize',15)
hold on

%% #2
load('SICSc05Eps0betaFAlpha80it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2)); 
cimin = zeros(1, sz(2));
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];
   
    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,6)
if isempty(g) == 0
    errorbar(bF(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 60])
xticks([0 20 40 60])
ylabel({'                           Mean time to pathogen extinction, \tau',''}, 'FontSize', 18)
set(gca,'Position',[0.1 0.1 0.14 0.27])
str = 'A(ii)';
dim = [.11 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on

%% #3
load('SICSc0Eps05betaFAlpha40it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2)); 
cimin = zeros(1, sz(2));
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];
   
    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,1)
if isempty(g) == 0
    errorbar(bF(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])
xlim([0 60])
xticks([0 20 40 60])
xticklabels([])
hold on

%% #4
load('SICSc05Eps0betaFAlpha40it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2)); 
cimin = zeros(1, sz(2));
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];
   
    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,1)
if isempty(g) == 0
    errorbar(bF(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])
xlim([0 60])
xticks([0 20 40 60])
xticklabels([])
set(gca,'Position',[0.1 0.41 0.14 0.27])
str = 'A(i)';
dim = [.11 .64 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on

%% #5                            Beta_D
load('SICSc0Eps05betaDAlpha80it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2));  
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak

    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];    

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,7)
if isempty(g) == 0
    errorbar(bD(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0.1 0.3])
xticks([0.1 0.2 0.3])
xlabel('\beta_D','FontSize',15)
hold on

%% #6
load('SICSc05Eps0betaDAlpha80it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2));  
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak

    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];    

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,7)
if isempty(g) == 0
    errorbar(bD(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0.1 0.3])
xticks([0.1 0.2 0.3])
set(gca,'Position',[0.28 0.1 0.14 0.27])
str = 'B(ii)';
dim = [.29 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on

%% #7
load('SICSc0Eps05betaDAlpha40it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2));  
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak

    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];    

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,2)
if isempty(g) == 0
    errorbar(bD(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])
xlim([0.1 0.3])
xticks([0.1 0.2 0.3])
xticklabels([])
hold on

%% #8
load('SICSc05Eps0betaDAlpha40it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2));  
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak

    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];    

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,2)
if isempty(g) == 0
    errorbar(bD(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])
xlim([0.1 0.3])
xticks([0.1 0.2 0.3])
xticklabels([])
set(gca,'Position',[0.28 0.41 0.14 0.27])
str = 'B(i)';
dim = [.29 .64 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on

%% #9                            Gamma
load('SICSc0Eps05gammaAlpha80it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2)); 
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,8)
if isempty(g) == 0
    errorbar(gam(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 80])
xticks([0 20 40 60 80])
xlabel('\gamma','FontSize',15)
hold on

%% #10
load('SICSc05Eps0gammaAlpha80it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2)); 
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,8)
if isempty(g) == 0
    errorbar(gam(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 80])
xticks([0 20 40 60 80])
set(gca,'Position',[0.46 0.1 0.14 0.27])
str = 'C(ii)';
dim = [.47 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on

%% #11
load('SICSc0Eps05gammaAlpha40it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2)); 
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,3)
if isempty(g) == 0
    errorbar(gam(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 50])
xlim([0 80])
xticks([0 20 40 60 80])
xticklabels([])
hold on

%% #12
load('SICSc05Eps0gammaAlpha40it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2)); 
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,3)
if isempty(g) == 0
    errorbar(gam(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 50])
xlim([0 80])
xticks([0 20 40 60 80])
xticklabels([])
set(gca,'Position',[0.46 0.41 0.14 0.27])
str = 'C(i)';
dim = [.56 .64 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on

%% #13                           Kappa
load('SICSc0Eps05kappaAlpha80it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2)); 
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,9)
if isempty(g) == 0
    errorbar(kap(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 160])
xticks([0 40 80 120 160])
xlabel('\kappa','FontSize',15)
hold on

%% #14
load('SICSc05Eps0kappaAlpha80it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2)); 
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,9)
if isempty(g) == 0
    errorbar(kap(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 160])
xticks([0 40 80 120 160])
set(gca,'Position',[0.64 0.1 0.14 0.27])
str = 'D(ii)';
dim = [.65 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on

%% #15
load('SICSc0Eps05kappaAlpha40it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2)); 
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,4)
if isempty(g) == 0
    errorbar(kap(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])
xlim([0 160])
xticks([0 40 80 120 160])
xticklabels([])
hold on

%% #16
load('SICSc05Eps0kappaAlpha40it100t500.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2)); 
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,4)
if isempty(g) == 0
    errorbar(kap(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])
xlim([0 160])
xticks([0 40 80 120 160])
xticklabels([])
set(gca,'Position',[0.64 0.41 0.14 0.27])
str = 'D(i)';
dim = [.65 .64 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on

%% #17                              Epsilon
T=1000;
load('SICSc05epsilonAlpha80it100t1000.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2)); 
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,10)
if isempty(g) == 0
    errorbar(eps(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 1])
xticks([0 0.5 1])
xlabel('\epsilon','FontSize',15)
set(gca,'Position',[0.82 0.1 0.14 0.27])
str = 'E(ii)';
dim = [.83 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on


%% #18
load('SICSc05epsilonAlpha40it100t1000.mat')
sz = size(minorextinct);
avgmin = zeros(1, sz(2));
cimin = zeros(1, sz(2)); 
for i = 1:sz(2)
    nzmin = nonzeros(minorextinct(:,i));  %all nonzero extinction times from a minor outbreak
    
    up2 = find(nzmin > T); %ignoring all extinction times greater than T
    nzmin(up2) = [];

    if length(nzmin) > 0.5*sz(1)  %finding the mean and CI for minor outbreaks
        avgmin(i) = mean(nzmin);
        cimin(i) = z*std(nzmin)/sqrt(length(nzmin));
    end
end
g = find(avgmin>0);

subplot(2,5,5)
if isempty(g) == 0
    errorbar(eps(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 120])
yticks([0 30 60 90 120])
xlim([0 1])
xticks([0 0.5 1])
xticklabels([])
set(gca,'Position',[0.82 0.41 0.14 0.27])
str = 'E(i)';
dim = [.83 .64 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on

