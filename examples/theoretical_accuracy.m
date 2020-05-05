%% Evaluate the PR/temporal theoretical accuracy of BSD.

P = struct; % Structure of generative model parameters.
P.a = 1; % Spike amplitude
P.sigma = 0.2; % Noise level
P.tauRise = 0.1; % Rise time
P.tauDecay = 0.5; % Decay time

O = struct; % Structure of options
O.dt = 0.1; % interval duration. (s)
O.superResolution = 1;
O.delta_max = 10; % limit lag: O.dt * O.delta_max; optional.
O.z1 = 2.366; % z1 and z2 parameters; optional.
O.z2 = 2.366;
O.Nsim = 1e4; % number of Monte Carlo runs; optional.
O.discretisation = 40; % discretisation parameter; optional.


[RF,DeltaT,lag,width] = BSD_theoretical_accuracy(P,O);

figure(1); hold on; plot(DeltaT,RF,'LineWidth',2); plot([0 0],[0,max(RF)],'LineWidth',2);
xlabel('Onset to spike (s)'); ylabel('Response function'); 
title(sprintf('Response function: delay to spike onset: %.3f s, width: %.3f s',lag,width));
set(gca,'FontSize',16);

% Notice a non-zero offset due to discretization error.

%% Show several noise levels in same plot
sigmas = [0,0.2,0.4];
for i = 1:length(sigmas);
    P.sigma = sigmas(i);
    [RF(:,i),DeltaT,lag,width] = BSD_theoretical_accuracy(P,O);
end;

figure(1); hold on; plot(DeltaT,RF,'LineWidth',2);
plot([0 0],[0,max(RF(:))],'LineWidth',2);
xlabel('Onset to spike (s)'); ylabel('Response function');  
legend(sprintf('\\sigma = %.2f',sigmas(1)),sprintf('\\sigma = %.2f',sigmas(2)),sprintf('\\sigma = %.2f',sigmas(3)))
title('Response functions');
set(gca,'FontSize',16);

%% Super-resolution.

P.sigma = 0.2;

O.superResolution = 1;
[RF1,DeltaT1] = BSD_theoretical_accuracy(P,O);
O.superResolution = 2;
[RF2,DeltaT2] = BSD_theoretical_accuracy(P,O);
O.superResolution = 3;
[RF3,DeltaT3] = BSD_theoretical_accuracy(P,O);
O.superResolution = 4;
[RF4,DeltaT4] = BSD_theoretical_accuracy(P,O);
O.superResolution = 5;
[RF5,DeltaT5] = BSD_theoretical_accuracy(P,O);


figure; hold on;
plot(DeltaT1,RF1,'LineWidth',2);
plot(DeltaT2,RF2,'LineWidth',2);
plot(DeltaT3,RF3,'LineWidth',2);
plot(DeltaT4,RF4,'LineWidth',2);
plot(DeltaT5,RF5,'LineWidth',2);
plot([0 0],[0,max(RF5(:))],'LineWidth',2);
xlabel('Onset to spike (s)'); ylabel('Response function');  
legend('No SR','SR=2','SR=3','SR=4','SR=5');
set(gca,'FontSize',16);

% Notice that the offset vanishes, and the width decreases.

%% Precision-Recall: True positive rate at fixed False Positive Rate


P = struct; % Structure of generative model parameters.
P.a = 1; % Spike amplitude
P.tauRise = 0.1; % Rise time
P.tauDecay = 0.5; % Decay time
nu = 0.1; % Spike Frequency.
O = struct; % Structure of options
O.superResolution = 1;
O.delta_max = 10; % limit lag: O.dt * O.delta_max; optional.
O.z1 = 2.366; % Usual z1
O.z2 = -2.366; % Minus sign here; means that we ensure \lambda_BSD = \lambda_1 (small false positive rate).
O.Nsim = 1e4; % number of Monte Carlo runs; increase for more accuracy
O.discretisation = 40; % discretisation parameter; optional.

fs = [1,2,3,4,5];
sigmas = [0:0.05:0.5];

TPR = zeros(length(fs),length(sigmas));
FPR = zeros(length(fs),length(sigmas));
FNR = zeros(length(fs),length(sigmas));

for i = 1:length(fs);
    display(fs(i));
    for j = 1:length(sigmas)
        O.dt = 1/fs(i);
        P.sigma = sigmas(j);
        [~,~,~,~,FPR(i,j),FNR(i,j)] = BSD_theoretical_accuracy(P,O);
        TPR(i,j) = 1- FNR(i,j);
        FPR(i,j) = FPR(i,j)/nu; % FPR is in false spikes/seconds. here, normalize to false spikes/ true spikes.
    end;
end;

Legend = cell(length(fs),1);
for i = 1:length(fs);
    Legend{i} = sprintf('f = %.f Hz',fs(i));
end;


figure;
plot(sigmas,TPR,'LineWidth',2);
xlabel('$\frac{\sigma}{a}$','interpreter','latex');
ylabel('True Positive Rate');
xlim([0,0.5]);
ylim([min(TPR(:)),1.02]);
legend(Legend,'location','SouthWest');
title('Spike Detection Probability, \tau_r = 0.1, \tau_d = 0.5')
set(gca,'FontSize',20);


%% Precision-Recall: best TPR/FNR trade-off using usual z1 and z2.


P = struct; % Structure of generative model parameters.
P.a = 1; % Spike amplitude
P.tauRise = 0.1; % Rise time
P.tauDecay = 0.5; % Decay time
nu = 0.1; % Spike Frequency.
O = struct; % Structure of options
O.superResolution = 1;
O.delta_max = 10; % limit lag: O.dt * O.delta_max; optional.
O.z1 = 2.366; % Usual z1
O.z2 = 2.366; % Usual z2
O.Nsim = 1e4; % number of Monte Carlo runs; increase for more accuracy
O.discretisation = 40; % discretisation parameter; optional.

fs = [1,2,3,4,5];
sigmas = [0:0.05:0.5];

TPR = zeros(length(fs),length(sigmas));
FPR = zeros(length(fs),length(sigmas));
FNR = zeros(length(fs),length(sigmas));

for i = 1:length(fs);
    display(fs(i));
    for j = 1:length(sigmas)
        O.dt = 1/fs(i);
        P.sigma = sigmas(j);
        [~,~,~,~,FPR(i,j),FNR(i,j)] = BSD_theoretical_accuracy(P,O);
        TPR(i,j) = 1- FNR(i,j);
        FPR(i,j) = FPR(i,j)/nu; % FPR is in false spikes/seconds. here, normalize to false spikes/ true spikes.
    end;
end;

Legend = cell(length(fs),1);
for i = 1:length(fs);
    Legend{i} = sprintf('f = %.f Hz',fs(i));
end;


figure('Position',[300,300,800,800]);
subplot(2,1,1);
plot(sigmas,TPR,'LineWidth',2); hold on;
xlabel('$\frac{\sigma}{a}$','interpreter','latex');
ylabel('True Positive Rate');
xlim([0,0.5]);
ylim([min(TPR(:)),1.02]);
set(gca,'FontSize',20);
title('Spike Detection Probability, \tau_r = 0.1, \tau_d = 0.5')
legend(Legend,'location','SouthWest');
subplot(2,1,2);
plot(sigmas,FPR,'LineWidth',2);
xlabel('$\frac{\sigma}{a}$','interpreter','latex');
ylabel('False Positive Rate');
xlim([0,0.5]);
set(gca,'FontSize',20);
