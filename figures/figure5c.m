%% Script for figure 5c.

addpath(genpath('../'))
load('../SpikeFinder/ground_truth_time_constants.mat');


datasets = [1,6,7,8,9,10];

tauRises = tauRiseGroundTruth(datasets);
tauDecays = tauDecayGroundTruth(datasets);

names = {'OGB-1','GCaMP5k','GCaMP6f','GCaMP6s','jRCAMP1a','jRGECO1a'};




fs = 10:5:200;
dts = 1./fs;
fref = 10;
sigma0s = [0,0.1,0.2,0.3];
sigmas = sqrt(fs/fref);

width = zeros(length(tauRises),length(sigmas),length(fs));
for i = 1:length(tauRises);
    for l = 1:length(sigma0s);
        for k = 1:length(fs);
            O = struct;
            O.dt = dts(k);
            O.z1 = 2.326;
            O.z2 = 2.326;
            O.delta_max = round(0.3/O.dt);
            O.discretisation = 40;        
            O.Nsim =1e5;
            P = struct;
            P.sigma = sigmas(k) * sigma0s(l);
            P.tauRise = tauRises(i);
            P.tauDecay = tauDecays(i);
            P.a = 1;
            P.b = 0;
        [~,~,~,width(i,k,l),~,~] = BSD_theoretical_accuracy(P,O);
        end;
    end;
end

save('data_figure5c.mat');

%%

load('data_figure5c.mat');

figure;
all_default_colors = get(gca,'colororder');
close;

colors = [4,2,1,3,6,7];


figure('Position',[200,200,1000,1200]);  hold on;
for i = 1: length(tauRises);
    subplot(2,3,i);
    plot(fs,fs,'LineWidth',2,'Color','Black');
    set(gca,'ColorOrderIndex',1)
    plot(fs,squeeze(1./width(i,:,:)),'LineWidth',2);    
    xlim([10,200]); ylim([10,200]);
    if i ==4
        leg=legend({'$\frac{\sigma}{a}_{10Hz} = 0$','$\frac{\sigma}{a}_{10Hz} = 0.1$', '$\frac{\sigma}{a}_{10Hz}  = 0.2$','$\frac{\sigma}{a}_{10Hz}  = 0.3$'},'location','northeast','Interpreter','Latex');
        legend boxoff
    end;
    set(gca,'FontSize',22);
    set(gcf, 'Color', 'w');
    xlabel('f (Hz)');
    ylabel('\delta t^{-1} (Hz)');
    title(sprintf('%s',names{i}));
end;
savefig('figure5c.fig');
export_fig('figure5c.png','-nocrop');
