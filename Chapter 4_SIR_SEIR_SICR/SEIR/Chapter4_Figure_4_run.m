%%% Producing Figure 4.4 in the thesis, chapter 4

%%% Pathogen extinction time plots with error bars, for varying model
%%% parameters

%%% This script is used in conjuction with the similarly named script found
%%% for SIR

%%% Run SIR script first followed by this


z = 1.96;    %95% confidence interval
bF = 0:10:60;
bD = (0:6)/20;
gam = 5:5:80;
kap = 10:10:160;
eps = 0:0.05:1;
T = 500;


%% Beta_F
load('SEIRbetaFAlpha80it100t500.mat')
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

subplot(3,5,11)
if isempty(g) == 0
    errorbar(bF(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])
xticks([0 20 40 60])
xlabel('\beta_F','FontSize',15)
set(gca,'Position',[0.1 0.1 0.14 0.27])
str = 'A(iii)';
dim = [.2 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on




load('SEIRbetaFAlpha40it100t500.mat')
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

subplot(3,5,6)
if isempty(g) == 0
    errorbar(bF(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])
xticks([0 20 40 60])
xticklabels([])
set(gca,'Position',[0.1 0.4 0.14 0.27])
ylabel('Mean time to pathogen extinction, \tau','FontSize',18)
str = 'A(ii)';
dim = [.2 .63 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on




load('SEIRbetaFAlpha10it100t500.mat')
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

subplot(3,5,1)
if isempty(g) == 0
    errorbar(bF(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])
xticks([0 20 40 60])
xticklabels([])
set(gca,'Position',[0.1 0.7 0.14 0.27])
str = 'A(i)';
dim = [.2 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on

%% Beta_D
load('SEIRbetaDAlpha80it100t500.mat')
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

subplot(3,5,12)
if isempty(g) == 0
    errorbar(bD(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xlim([0.05 0.3])
xticks([0.1 0.2 0.3])
xlabel('\beta_D','FontSize',15)
set(gca,'Position',[0.28 0.1 0.14 0.27])
str = 'B(iii)';
dim = [.37 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');


load('SEIRbetaDAlpha40it100t500.mat')
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

subplot(3,5,7)
if isempty(g) == 0
    errorbar(bD(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xlim([0.05 0.3])
xticks([0.1 0.2 0.3])
xticklabels([])
set(gca,'Position',[0.28 0.4 0.14 0.27])
str = 'B(ii)';
dim = [.37 .63 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');


load('SEIRbetaDAlpha1it100t500.mat')
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

subplot(3,5,2)
if isempty(g) == 0
    errorbar(bD(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xlim([0.05 0.3])
xticks([0.1 0.2 0.3])
xticklabels([])
set(gca,'Position',[0.28 0.7 0.14 0.27])
str = 'B(i)';
dim = [.37 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');

%% Gamma
load('SEIRgammaAlpha80it100t500.mat')
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

subplot(3,5,13)
if isempty(g) == 0
    errorbar(gam(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])
xlim([0 80])
xticks([0 20 40 60 80])
xlabel('\gamma','FontSize',15)
set(gca,'Position',[0.46 0.1 0.14 0.27])
str = 'C(iii)';
dim = [.56 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');


load('SEIRgammaAlpha40it100t500.mat')
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

subplot(3,5,8)
if isempty(g) == 0
    errorbar(gam(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])
xlim([0 80])
xticks([0 20 40 60 80])
xticklabels([])
set(gca,'Position',[0.46 0.4 0.14 0.27])
str = 'C(ii)';
dim = [.56 .63 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');



load('SEIRgammaAlpha10it100t500.mat')
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

subplot(3,5,3)
if isempty(g) == 0
    errorbar(gam(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])
xlim([0 80])
xticks([0 20 40 60 80])
xticklabels([])
set(gca,'Position',[0.46 0.7 0.14 0.27])
str = 'C(i)';
dim = [.56 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');

%% Kappa
load('SEIRkappaAlpha80it100t500.mat')
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

subplot(3,5,14)
if isempty(g) == 0
    errorbar(kap(g), avgmin(g), cimin(g), 'g-')
    hold on
end
ylim([0 20])
xlim([0 160])
xticks([0 40 80 120 160])
set(gca,'Position',[0.64 0.1 0.14 0.27])
xlabel('\kappa','FontSize',15)
str = 'D(iii)';
dim = [.74 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');


load('SEIRkappaAlpha40it100t500.mat')
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

subplot(3,5,9)
if isempty(g) == 0
    errorbar(kap(g), avgmin(g), cimin(g), 'g-')
    hold on
end
ylim([0 20])
xlim([0 160])
xticks([0 40 80 120 160])
set(gca,'Position',[0.64 0.4 0.14 0.27])
xticklabels([])
str = 'D(ii)';
dim = [.74 .63 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');


load('SEIRkappaAlpha10it100t500.mat')
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

subplot(3,5,4)
if isempty(g) == 0
    errorbar(kap(g), avgmin(g), cimin(g), 'g-')
    hold on
end
ylim([0 20])
xlim([0 160])
xticks([0 40 80 120 160])
set(gca,'Position',[0.64 0.7 0.14 0.27])
xticklabels([])
str = 'D(i)';
dim = [.74 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');

%% Epsilon
load('SEIRepsilonAlpha80it100t500.mat')
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

subplot(3,5,15)
if isempty(g) == 0
    errorbar(eps(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
set(gca,'Position',[0.82 0.1 0.14 0.27])
xlabel('\epsilon','FontSize',15)
str = 'E(iii)';
dim = [.92 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');



load('SEIRepsilonAlpha40it100t500.mat')
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

subplot(3,5,10)
if isempty(g) == 0
    errorbar(eps(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
set(gca,'Position',[0.82 0.4 0.14 0.27])
xticklabels([])
str = 'E(ii)';
dim = [.92 .63 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');



load('SEIRepsilonAlpha10it100t500.mat')
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

subplot(3,5,5)
if isempty(g) == 0
    errorbar(eps(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xticklabels([])
set(gca,'Position',[0.82 0.7 0.14 0.27])
str = 'E(i)';
dim = [.92 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');