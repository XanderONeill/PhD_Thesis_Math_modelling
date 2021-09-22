%  Chapter 2: Plotting the relevant regions of transmission coefficients
%  satisfying the epidemiological criteria for three different values of
%  rho. Relevant Figure: Fig A.2

%   Uses matrices 'varallrho083', 'varallrho085' and 'varallrho090' which
%   were acquired by running the code titled

%% 2D plot when rho = 0.83
clearvars
load('varallrho083')
var = varintersect;
a = unique(var(:,1),'sorted');
upper = [];
lower = [];
for i = 1:length(a)
    ind = find(var(:,1) == a(i));
    upper = [upper, max(var(ind(1):ind(end),2))];
    lower = [lower, min(var(ind(1):ind(end),2))];   
end

figure
subplot(3,1,1)
plot(a, upper, 'k-', a, lower, 'k-');
hold on
b = [a, a(end:-1:1)];
b = b(:);
inBetween = [lower, fliplr(upper)];
fill(b, inBetween, 'k');
title('A', 'Units', 'normalized', 'Position',[0.95 0.82])
ylabel('\beta_E');
ylim([1.15 2.35])
xlim([61.5 67.1])
set(gca, 'Xticklabels', []);
set(gca, 'Position', [0.08 0.7 0.85, 0.27]);
hold off


%% 2D Plot when rho = 0.85

clearvars
load('varallrho085')
var = varintersect;
a = unique(var(:,1),'sorted');
upper = [];
lower = [];
for i = 1:length(a)
    ind = find(var(:,1) == a(i));
    upper = [upper, max(var(ind(1):ind(end),2))];
    lower = [lower, min(var(ind(1):ind(end),2))];   
end

subplot(3,1,2)
plot(a, upper, 'k-', a, lower, 'k-');
hold on
b = [a, a(end:-1:1)];
b = b(:);
inBetween = [lower, fliplr(upper)];
fill(b, inBetween, 'k');
title('B', 'Units', 'normalized', 'Position',[0.95 0.82])
ylabel('\beta_E');
ylim([1.15 2.35])
xlim([61.5 67.1])
set(gca, 'Xticklabels', []);
set(gca, 'Position', [0.08 0.4 0.85, 0.27]);
hold off

%% 2D Plot when rho = 0.90

clearvars
load('varallrho090')
var = varintersect;
a = unique(var(:,1),'sorted');
upper = [];
lower = [];
for i = 1:length(a)
    ind = find(var(:,1) == a(i));
    upper = [upper, max(var(ind(1):ind(end),2))];
    lower = [lower, min(var(ind(1):ind(end),2))];   
end

subplot(3,1,3)
plot(a, upper, 'k-', a, lower, 'k-');
hold on
b = [a, a(end:-1:1)];
b = b(:);
inBetween = [lower, fliplr(upper)];
fill(b, inBetween, 'k');
title('C', 'Units', 'normalized', 'Position',[0.95 0.82])
ylabel('\beta_E');
xlabel('\beta_F')
ylim([1.15 2.35])
xlim([61.5 67.1])
set(gca, 'Position', [0.08 0.1 0.85, 0.27]);
hold off
