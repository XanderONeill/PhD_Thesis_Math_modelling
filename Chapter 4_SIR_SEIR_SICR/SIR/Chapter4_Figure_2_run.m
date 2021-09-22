%%% Producing Figure 4.2 in the thesis, chapter 4
%%% Pathogen extinction time plots with error bars, for a varying disease 
%%% induced mortality rate \alpha

%%% This script is used in conjuction with the similarly named script found
%%% for SEIS, SIR and SEIR

%%% Run in order of SIS, SEIS, SIR, SEIR

z = 1.96;    %95% confidence interval
alph = 5:5:80;
T = 500;

%% SIRAlpha
load('SIRAlphait100t500.mat')

sz = size(minorextinct);
avgmin = zeros(1, sz(2));  cimin = zeros(1, sz(2)); 
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

subplot(1,2,2)
if isempty(g) == 0
    errorbar(alph(g), avgmin(g), cimin(g), 'k-', 'LineWidth', 0.5)
    hold on
end
xlim([0 80])
