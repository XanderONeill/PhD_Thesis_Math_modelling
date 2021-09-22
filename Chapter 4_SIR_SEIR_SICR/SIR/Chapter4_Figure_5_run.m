%%% Producing Figure 4.5 in the thesis, chapter 4
%%% Pathogen extinction time plots. No error bars needed as the varied
%%% parameter has no impact on the SIR model. This script is just used to
%%% produce the straight black lines seen in Fig 4.5

%%% This script is used in conjuction with the similarly named script found
%%% for SICS, SIR and SICR

%%% Run in order of SIS, SICS, SIR, SICR

T = 500;
cvec = 0:0.05:1;

load('SIRbetaDAlpha80it100t500')
nz = nonzeros(minorextinct(:,4));
up = nz > T; %ignoring all extinction times greater than T
nz(up) = [];
avg = mean(nz);

avgvec = ones(length(cvec),1)*avg;
subplot(2,2,4)
plot(cvec, avgvec,'k-','LineWidth',1)
hold on


load('SIRbetaDAlpha40it100t500')
nz = nonzeros(minorextinct(:,4));
up = find(nz > T); %ignoring all extinction times greater than T
nz(up) = [];
avg = mean(nz);

avgvec = ones(length(cvec),1)*avg;
subplot(2,2,2)
plot(cvec, avgvec,'k-','LineWidth',1)
hold on