%%% Script to produce Figure 5.7 in the thesis, Chapter 5

%%% 2x3 cross-sectional plot for a fixed small mammal density and varying
%%% ungulate density, or a fixed ungulate density and varying small mammal
%%% density. This is done for a low, medium and high fixed small mammal
%%% density along one row, and low, medium and high fixed ungulate density 
%%% along the other.

%%% uses the both the functions (models):
%%%     Model_Ext_Rosa_Pugliese and Model_tick_host
%%% and compares the cross sectional tick densities for each

clearvars
close all

%% Running the model for fixed ungulate densities

Hsvec = [1:10,20:5:500]; %varying ungulate densities
Huvec = [10, 20, 40];      %fixed low, medium and high ungulate densities
LqR = zeros(length(Hsvec),3); LfR = zeros(length(Hsvec),3); %3 columns for low, medium and high ungulate levels
NqR = zeros(length(Hsvec),3); NfR = zeros(length(Hsvec),3); 
AqR = zeros(length(Hsvec),3); AfR = zeros(length(Hsvec),3); 

LqM = zeros(length(Hsvec),3); LfM = zeros(length(Hsvec),3); 
NqM = zeros(length(Hsvec),3); NfsM = zeros(length(Hsvec),3); NfuM = zeros(length(Hsvec),3);
AqM = zeros(length(Hsvec),3); AfM = zeros(length(Hsvec),3); 

for j = 1:3
    Hu = Huvec(j);
    for i = 1:length(Hsvec)
        Hs = Hsvec(i);

        y0 = 100*ones(6,1);
        options = odeset('Refine', 1,'RelTol',1e-5,'AbsTol',1e-5,'NonNegative',1:6);
        [t,y] = ode15s(@(t,y) Model_Ext_Rosa_Pugliese(t, y, Hs, Hu), [0 10000], y0, options);
        LqR(i,j) = y(end,1); LfR(i,j) = y(end,2); NqR(i,j) = y(end,3); 
        NfR(i,j) = y(end,4); AqR(i,j) = y(end,5); AfR(i,j) = y(end,6);

        y0 = 100*ones(7,1);
        options = odeset('Refine', 1,'RelTol',1e-5,'AbsTol',1e-5,'NonNegative',1:7);
        [t1,y1] = ode15s(@(t1,y1) Model_tick_host(t1, y1, Hs, Hu), [0 10000], y0, options);
        LqM(i,j) = y1(end,1); LfM(i,j) = y1(end,2); NqM(i,j) = y1(end,3); 
        NfsM(i,j) = y1(end,4); NfuM(i,j) = y1(end,5); AqM(i,j) = y1(end,6); 
        AfM(i,j) = y1(end,7);
    end
end

%%% Plotting the cross-sectionals
figure
subplot(3,2,1)
set(gca, 'Position', [0.08 0.68 0.27 0.25]);
title('A(i): H_u = 10')
set(gca, 'XTickLabels',[]);
xticks([0 250 500]);
xlim([0 500]);
yyaxis left
plot(Hsvec, LfR(:,1), Hsvec, LfM(:,1), ':', 'LineWidth', 1.5)
yyaxis right
plot(Hsvec, AfR(:,1), Hsvec, AfM(:,1), ':','LineWidth', 1.5)

subplot(3,2,3)
set(gca, 'Position', [0.08 0.38 0.27 0.25]);
title('A(ii): H_u = 20')
set(gca, 'XTickLabels',[]);
xticks([0 250 500]);
xlim([0 500]);
yyaxis left
plot(Hsvec, LfR(:,2), Hsvec, LfM(:,2), ':','LineWidth', 1.5)
ylabel('Feeding Larvae Density')
yyaxis right
plot(Hsvec, AfR(:,2), Hsvec, AfM(:,2), ':','LineWidth', 1.5)

subplot(3,2,5)
set(gca, 'Position', [0.08 0.08 0.27 0.25]);
title('A(iii): H_u = 40');
xlabel('Small Mammal Density')
xticks([0 250 500]);
xlim([0 500]);
yyaxis left
plot(Hsvec, LfR(:,3), Hsvec, LfM(:,3), ':','LineWidth', 1.5)
yyaxis right
plot(Hsvec, AfR(:,3), Hsvec, AfM(:,3), ':','LineWidth', 1.5)

%% Running the model for fixed small mammal densities

Huvec = [1:10,20:5:250];
Hsvec = [42, 85, 170];
LqRu = zeros(length(Huvec),3); LfRu = zeros(length(Huvec),3); 
NqRu = zeros(length(Huvec),3); NfRu = zeros(length(Huvec),3); 
AqRu = zeros(length(Huvec),3); AfRu = zeros(length(Huvec),3); 

LqMu = zeros(length(Huvec),3); LfMu = zeros(length(Huvec),3); 
NqMu = zeros(length(Huvec),3); NfsMu = zeros(length(Huvec),3); NfuMu = zeros(length(Huvec),3);
AqMu = zeros(length(Huvec),3); AfMu = zeros(length(Huvec),3); 

for j = 1:3
    Hs = Hsvec(j);
    for i = 1:length(Huvec)
        Hu = Huvec(i);

        y0 = 100*ones(6,1);
        options = odeset('Refine', 1,'RelTol',1e-5,'AbsTol',1e-5,'NonNegative',1:6);
        [t,y] = ode15s(@(t,y) Model_Ext_Rosa_Pugliese(t, y, Hs, Hu), [0 10000], y0, options);
        LqRu(i,j) = y(end,1); LfRu(i,j) = y(end,2); NqRu(i,j) = y(end,3); 
        NfRu(i,j) = y(end,4); AqRu(i,j) = y(end,5); AfRu(i,j) = y(end,6);

        y0 = 100*ones(7,1);
        options = odeset('Refine', 1,'RelTol',1e-5,'AbsTol',1e-5,'NonNegative',1:7);
        [t1,y1] = ode15s(@(t1,y1) Model_tick_host(t1, y1, Hs, Hu), [0 10000], y0, options);
        LqMu(i,j) = y1(end,1); LfMu(i,j) = y1(end,2); NqMu(i,j) = y1(end,3); 
        NfsMu(i,j) = y1(end,4); NfuMu(i,j) = y1(end,5); AqMu(i,j) = y1(end,6); 
        AfMu(i,j) = y1(end,7);
    end
end

%%% Producing the respective cross-sectional plots
subplot(3,2,2)
set(gca, 'Position', [0.45 0.68 0.27 0.25]);
title('B(i): H_s = 42');
set(gca, 'XTickLabels',[]);
xticks([0 125 250]);
xlim([0 250]);
yyaxis left
plot(Huvec, LfRu(:,1), Huvec, LfMu(:,1), ':','LineWidth', 1.5)
yyaxis right
plot(Huvec, AfRu(:,1), Huvec, AfMu(:,1), ':','LineWidth', 1.5)

subplot(3,2,4)
set(gca, 'Position', [0.45 0.38 0.27 0.25]);
title('B(ii): H_s = 85');
set(gca, 'XTickLabels',[]);
xticks([0 125 250]);
xlim([0 250]);
yyaxis left
plot(Huvec, LfRu(:,2), Huvec, LfMu(:,2), ':','LineWidth', 1.5)
yyaxis right
plot(Huvec, AfRu(:,2), Huvec, AfMu(:,2), ':','LineWidth', 1.5)
ylabel({'Feeding Adult Density'});

subplot(3,2,6)
set(gca, 'Position', [0.45 0.08 0.27 0.25]);
title('B(iii): H_s = 170');
xlabel('Ungulate Density')
xticks([0 125 250]);
xlim([0 250]);
yyaxis left
plot(Huvec, LfRu(:,3), Huvec, LfMu(:,3), ':','LineWidth', 1.5)
yyaxis right
plot(Huvec, AfRu(:,3), Huvec, AfMu(:,3), ':','LineWidth', 1.5)