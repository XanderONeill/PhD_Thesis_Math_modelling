%%% Producing Figure 4.3 in the thesis, chapter 4

%%% Pathogen extinction time plots with error bars, for varying model
%%% parameters

%%% This script is used in conjuction with the similarly named script found
%%% for SIS

%%% Run the SIS script first followed by this


z = 1.96;    %95% confidence interval
bF = 0:10:60;
bD = (0:6)/20;
gam = 5:5:80;
kap = 10:10:160;
eps = 0:0.05:1;
T = 500;


%% Beta_F
load('SEISbetaFAlpha80it100t500.mat')
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
    errorbar(bF(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
set(gca,'Position',[0.1 0.41 0.14 0.27])
xticks([0 20 40 60])
xticklabels([])
%ylim([0 2])
hold on
str = 'A(i)';
dim = [.2 .64 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');

%% Beta_D
load('SEISbetaDAlpha80it100t500.mat')
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
    errorbar(bD(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
set(gca,'Position',[0.28 0.41 0.14 0.27])
xticks([0.15 0.2 0.25 0.3])
xticklabels([])
str = 'B(i)';
dim = [.38 .64 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
%ylim([0 2])

%% Gamma
load('SEISgammaAlpha80it100t500.mat')
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
    errorbar(gam(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
set(gca,'Position',[0.46 0.41 0.14 0.27])
xticks([0 40 80])
xticklabels([])
str = 'C(i)';
dim = [.56 .64 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
%ylim([0 2])

%% Kappa
load('SEISkappaAlpha80it100t500.mat')
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
    errorbar(kap(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
set(gca,'Position',[0.64 0.41 0.14 0.27])
xticks([0 40 80 120 160])
xlim([0 160])
xticklabels([])
str = 'D(i)';
dim = [.74 .64 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
%ylim([0 2])

%% Epsilon
T = 1000;
load('SEISepsilonAlpha80it100t1000.mat')
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
    errorbar(eps(g), avgmin(g), cimin(g), 'g-', 'LineWidth', 0.5)
    hold on
end
set(gca,'Position',[0.82 0.41 0.14 0.27])
xticklabels([])
str = 'E(i)';
dim = [.92 .64 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
%ylim([0 2])
