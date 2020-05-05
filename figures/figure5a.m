%% Script for Figure 5a: TPR vs sigma.


fs = [1,2,3,4,5];
sigmas = [0:0.01:0.5];

a = 1;
tauRise = 0.1;
tauDecay = 0.5;
nu = 0.1;


TPR = zeros(length(fs),length(sigmas));
FDR = zeros(length(fs),length(sigmas));
tilde_sigma = zeros(length(fs),length(sigmas));

for i = 1:length(fs);
    for j = 1:length(sigmas)
        dt = 1/fs(i);
        sigma = sigmas(j);
        tilde_sigma(i,j) = sigma/sqrt(dt);
        P = struct;
        P.tauRise = tauRise;
        P.tauDecay  = tauDecay;
        P.a = 1;
        P.sigma = sigma;
        O = struct; O.dt = dt; O.nNeurons=1;
        O.z1 = 2.366;
        O.z2 = -z1; % Make sure that lambda_BSD = lambda_1, i.e. set a fixed False Positive Rate.
        O.Nsim = 1e5; % number of Monte Carlo runs; optional., default = 1e4
        [~,~,~,~,FPR,FNR] = BSD_theoretical_accuracy(P,O);
        TPR(i,j)= 1-FNR;
        FDR(i,j) = FPR/nu; % FPR is in false spikes/seconds. here, normalize to false spikes/ true spikes.
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
savefig('figure5a.fig');
set(gcf,'Color','White');
export_fig('figure5a.png','-nocrop');


