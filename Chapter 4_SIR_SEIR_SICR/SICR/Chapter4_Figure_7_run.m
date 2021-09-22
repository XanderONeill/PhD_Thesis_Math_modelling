%%% Producing Figure 4.7 in the thesis, chapter 4

%%% Pathogen extinction time plots with error bars, for varying model
%%% parameters

%%% This script is used in conjuction with the similarly named script found
%%% for SIR

%%% Run the SIR script followed by this

z = 1.96;    %95% confidence interval
bF = 0:10:60;
bD = (0:6)/20;
gam = 5:5:80;
kap = 10:10:160;
eps = 0:0.05:1;
T = 500;

%%%plots 1-6 varying betaF: ODD (c=0,epsilon=0.5) EVEN (c=0.5,epsilon=0.5)
%%% 7-12 varying betaD: ODD (c=0,epsilon=0.5) EVEN (c=0.5,epsilon=0.5)
%%% 13-18 varying gamma: ODD (c=0,epsilon=0.5) EVEN (c=0.5,epsilon=0.5)
%%% 19-24 varying kappa: ODD (c=0,epsilon=0.5) EVEN (c=0.5,epsilon=0.5)
%%% 25-27 varying epsilon: C=0.5

%% #1:                      beta_F

load('SICRc0Eps05betaFAlpha80it100t500.mat')
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
    errorbar(bF(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xticks([0 20 40 60])
xlabel('\beta_F', 'FontSize',15)
hold on

%% #2
load('SICRc05betaFAlpha80it100t500.mat')
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
    errorbar(bF(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xticks([0 20 40 60])
hold off



%% #3
load('SICRc0Eps05betaFAlpha40it100t500.mat')
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
    errorbar(bF(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xticks([0 20 40 60])
xticklabels([])
ylabel({'Mean time to pathogen extinction, \tau',''},'FontSize',18)
hold on

%% #4
load('SICRc05betaFAlpha40it100t500.mat')
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
    errorbar(bF(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xticks([0 20 40 60])
xticklabels([])
hold on


%% #5
load('SICRc0Eps05betaFAlpha1it100t500.mat')
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
    errorbar(bF(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xticks([0 20 40 60])
xticklabels([])
hold on

%% #6
load('SICRc05betaFAlpha1it100t500.mat')
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
    errorbar(bF(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xticks([0 20 40 60])
xticklabels([])
hold on

%% #7                           Beta_D
load('SICRc0Eps05betaDAlpha80it100t500.mat')
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
    errorbar(bD(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xlim([0.05 0.3])
xticks([0.1 0.2 0.3])
xlabel('\beta_D', 'FontSize',15)

%% #8 
load('SICRc05betaDAlpha80it100t500.mat')
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
    errorbar(bD(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xlim([0.05 0.3])
xticks([0.1 0.2 0.3])


%% #9 
load('SICRc0Eps05betaDAlpha40it100t500.mat')
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
    errorbar(bD(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xlim([0.05 0.3])
xticks([0.1 0.2 0.3])
xticklabels([])

%% #10 
load('SICRc05betaDAlpha40it100t500.mat')
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
    errorbar(bD(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xlim([0.05 0.3])
xticks([0.1 0.2 0.3])
xticklabels([])



%% #11 
load('SICRc0Eps05betaDAlpha1it100t500.mat')
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
    errorbar(bD(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xlim([0.05 0.3])
xticks([0.1 0.2 0.3])
xticklabels([])

%% #12 
load('SICRc05betaDAlpha1it100t500.mat')
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
    errorbar(bD(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 4])
xlim([0.05 0.3])
xticks([0.1 0.2 0.3])
xticklabels([])

%% #13                              Gamma
load('SICRc0Eps05gammaAlpha80it100t500.mat')
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
    errorbar(gam(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 80])
xticks([0 20 40 60 80])
xlabel('\gamma', 'FontSize',15)

%% #14
load('SICRc05gammaAlpha80it100t500.mat')
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
    errorbar(gam(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 80])
xticks([0 20 40 60 80])

%% #15
load('SICRc0Eps05gammaAlpha40it100t500.mat')
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
    errorbar(gam(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 80])
xticks([0 20 40 60 80])
xticklabels([])

%% #16
load('SICRc05gammaAlpha40it100t500.mat')
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
    errorbar(gam(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 80])
xticks([0 20 40 60 80])
xticklabels([])

%% #17
load('SICRc0Eps05gammaAlpha1it100t500.mat')
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
    errorbar(gam(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 80])
xticks([0 20 40 60 80])
xticklabels([])

%% #18
load('SICRc05gammaAlpha1it100t500.mat')
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
    errorbar(gam(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 80])
xticks([0 20 40 60 80])
xticklabels([])

%% #19                               Kappa
load('SICRc0Eps05kappaAlpha80it100t500.mat')
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
    errorbar(kap(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 160])
xticks([0 40 80 120 160])
xlabel('\kappa', 'FontSize',15)

%% #20
load('SICRc05kappaAlpha80it100t500.mat')
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
    errorbar(kap(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 160])
xticks([0 40 80 120 160])

%% #21
load('SICRc0Eps05kappaAlpha40it100t500.mat')
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
    errorbar(kap(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 160])
xticks([0 40 80 120 160])
xticklabels([])

%% #22
load('SICRc05kappaAlpha40it100t500.mat')
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
    errorbar(kap(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 160])
xticks([0 40 80 120 160])
xticklabels([])

%% #23
load('SICRc0Eps05kappaAlpha1it100t500.mat')
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
    errorbar(kap(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 160])
xticks([0 40 80 120 160])
xticklabels([])

%% #24
load('SICRc05kappaAlpha1it100t500.mat')
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
    errorbar(kap(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 5])
xlim([0 160])
xticks([0 40 80 120 160])
xticklabels([])

%% #25                               Epsilon
T=1000;
load('SICRc05epsilonAlpha80it100t1000.mat')
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
    errorbar(eps(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])
xlim([0 1])
xticks([0 0.5 1])
xlabel('\epsilon', 'FontSize',15)

%% #26
load('SICRc05epsilonAlpha40it100t1000.mat')
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
    errorbar(eps(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])
xlim([0 1])
xticks([0 0.5 1])
xticklabels([])

%% #27
load('SICRc05epsilonAlpha1it100t500.mat')
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
    errorbar(eps(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])
xlim([0 1])
xticks([0 0.5 1])
xticklabels([])


%% Tightening of Plots

subplot(3,5,1)
set(gca, 'Position',[0.1 0.7 0.14 0.27])
str = 'A(i)';
dim = [.2 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
subplot(3,5,6)
set(gca, 'Position',[0.1 0.4 0.14 0.27])
str = 'A(ii)';
dim = [.2 .63 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
subplot(3,5,11)
set(gca, 'Position',[0.1 0.1 0.14 0.27])
str = 'A(iii)';
dim = [.2 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');

subplot(3,5,2)
set(gca, 'Position',[0.28 0.7 0.14 0.27])
str = 'B(i)';
dim = [.38 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
subplot(3,5,7)
set(gca, 'Position',[0.28 0.4 0.14 0.27])
str = 'B(ii)';
dim = [.38 .63 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
subplot(3,5,12)
set(gca, 'Position',[0.28 0.1 0.14 0.27])
str = 'B(iii)';
dim = [.38 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');

subplot(3,5,3)
set(gca, 'Position',[0.46 0.7 0.14 0.27])
str = 'C(i)';
dim = [.56 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
subplot(3,5,8)
set(gca, 'Position',[0.46 0.4 0.14 0.27])
str = 'C(ii)';
dim = [.56 .63 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
subplot(3,5,13)
set(gca, 'Position',[0.46 0.1 0.14 0.27])
str = 'C(iii)';
dim = [.56 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');

subplot(3,5,4)
set(gca, 'Position',[0.64 0.7 0.14 0.27])
str = 'D(i)';
dim = [.74 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
subplot(3,5,9)
set(gca, 'Position',[0.64 0.4 0.14 0.27])
str = 'D(ii)';
dim = [.74 .63 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
subplot(3,5,14)
set(gca, 'Position',[0.64 0.1 0.14 0.27])
str = 'D(iii)';
dim = [.74 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');

subplot(3,5,5)
set(gca, 'Position',[0.82 0.7 0.14 0.27])
str = 'E(i)';
dim = [.92 .93 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
subplot(3,5,10)
set(gca, 'Position',[0.82 0.4 0.14 0.27])
str = 'E(ii)';
dim = [.92 .63 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
subplot(3,5,15)
set(gca, 'Position',[0.82 0.1 0.14 0.27])
str = 'E(iii)';
dim = [.92 .33 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');

