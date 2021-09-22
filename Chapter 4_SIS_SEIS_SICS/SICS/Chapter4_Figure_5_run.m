%%% Producing Figure 4.5 from the thesis, Chapter 4.

%%% Producing plots with error bars for a varying chronically infected 
%%% death rate alpha_C

%%% This is done for epsilon = 0, epsilon = 0.5. 

%%% Run in conjunction with the similarly named script in the SIS, SIR and 
%%% SICR folder.

%%% Run in order of SIS, SICS, SIR, SICR

z = 1.96;
T = 500;

cvec = 0:0.05:1;
%% Alpha = 80, Epsilon = 0
load('SICScEps0Alpha80it100t500')
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

subplot(2,2,3)
if isempty(g) == 0
    errorbar(cvec(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end
ylim([0 3])

%% Alpha = 40, Epsilon = 0
load('SICScEps0Alpha40it100t500')
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

subplot(2,2,1)
if isempty(g) == 0
    errorbar(cvec(g), avgmin(g), cimin(g), 'b-', 'LineWidth', 0.5)
    hold on
end

%% Alpha = 80, Epsilon 0.5
load('SICScEps05Alpha80it100t500')
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

subplot(2,2,3)
if isempty(g) == 0
    errorbar(cvec(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 3])
xlim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1])
xlabel('c');
ylabel({'                                               Mean time to pathogen extinction, \tau',''},'FontSize',15)
set(gca, 'Position', [0.27 0.1 0.21 0.405])
str = 'A(ii)';
dim = [.44 .465 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on

%% Alpha = 40, Epsilon = 0.5
load('SICScEps05Alpha40it100t500')
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

subplot(2,2,1)
if isempty(g) == 0
    errorbar(cvec(g), avgmin(g), cimin(g), 'r-', 'LineWidth', 0.5)
    hold on
end
ylim([0 100])
xlim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1])
xticklabels([])
set(gca, 'Position', [0.27 0.55 0.21 0.405])
str = 'A(i)';
dim = [.44 .915 .03 .03];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'off', 'EdgeColor', 'w');
hold on