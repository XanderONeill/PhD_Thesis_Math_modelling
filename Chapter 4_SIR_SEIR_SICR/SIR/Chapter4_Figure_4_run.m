%%% Producing Figure 4.4 in the thesis, chapter 4

%%% Pathogen extinction time plots with error bars, for varying model
%%% parameters

%%% This script is used in conjuction with the similarly named script found
%%% for SEIR

%%% Run this script first followed by the SEIR script

z = 1.96;    %95% confidence interval
bF = 0:10:60;
bD = (0:6)/20;
gam = 5:5:80;
T = 500;

figure

%% Beta_F
load('SIRbetaFAlpha80it100t500.mat')
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
    errorbar(bF(g), avgmin(g), cimin(g), 'k-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])
hold on


load('SIRbetaFAlpha40it100t500.mat')
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
    errorbar(bF(g), avgmin(g), cimin(g), 'k-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])
hold on



load('SIRbetaFAlpha10it100t500.mat')
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
    errorbar(bF(g), avgmin(g), cimin(g), 'k-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])
hold on

%% Beta_D
load('SIRbetaDAlpha80it100t500.mat')
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
    errorbar(bD(g), avgmin(g), cimin(g), 'k-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])



load('SIRbetaDAlpha40it100t500.mat')
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
    errorbar(bD(g), avgmin(g), cimin(g), 'k-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])



load('SIRbetaDAlpha10it100t500.mat')
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
    errorbar(bD(g), avgmin(g), cimin(g), 'k-', 'LineWidth', 0.5)
    hold on
end
ylim([0 2])
%% Gamma
load('SIRgammaAlpha80it100t500.mat')
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
    errorbar(gam(g), avgmin(g), cimin(g), 'k-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])

Kapvec80 = ones(16,1)*avgmin(10); %for later use in kappa and epsilon plots
Epsvec80 = ones(21,1)*avgmin(10);

load('SIRgammaAlpha40it100t500.mat')
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
    errorbar(gam(g), avgmin(g), cimin(g), 'k-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])

Kapvec40 = ones(16,1)*avgmin(10); %for later use in kappa and epsilon plots
Epsvec40 = ones(21,1)*avgmin(10);


load('SIRgammaAlpha10it100t500.mat')
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
    errorbar(gam(g), avgmin(g), cimin(g), 'k-', 'LineWidth', 0.5)
    hold on
end
ylim([0 20])


Kapvec10 = ones(16,1)*avgmin(10); %avgmin(10) corresponds to gamma = 50
Epsvec10 = ones(21,1)*avgmin(10); %which is baseline


subplot(3,5,4)
plot(10:10:160,Kapvec10,'k-');
hold on
subplot(3,5,9)
plot(10:10:160,Kapvec40,'k-');
hold on
subplot(3,5,14)
plot(10:10:160,Kapvec80,'k-');
hold on

subplot(3,5,5)
plot(0:0.05:1,Epsvec10,'k-');
hold on
subplot(3,5,10)
plot(0:0.05:1,Epsvec40,'k-');
hold on
subplot(3,5,15)
plot(0:0.05:1,Epsvec80,'k-');
hold on