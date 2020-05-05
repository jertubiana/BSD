load('ground_truth_time_constants.mat');

all_aInferred_groundTruth = all_aInferred;
all_tauRiseInferred_groundTruth = all_tauRiseInferred;
all_tauDecayInferred_groundTruth = all_tauDecayInferred;

%%
folder = '../SpikeFinder/predictions_blind/';


all_aInferred_blind = cell(10,1);
all_tauRiseInferred_blind = cell(10,1);
all_tauDecayInferred_blind = cell(10,1);


for dataset = 1:10
    load(strcat(folder,sprintf('%d.train.mat',dataset) ) )
    aInferred_blind = zeros(nNeurons,1);
    tauRiseInferred_blind = zeros(nNeurons,1);
    tauDecayInferred_blind = zeros(nNeurons,1);
    
    for i = 1:nNeurons;
        aInferred_blind(i) = pAlg{i}.a;
        tauRiseInferred_blind(i) = pAlg{i}.tauRise;
        tauDecayInferred_blind(i) = pAlg{i}.tauDecay;
    end;
    all_aInferred_blind{dataset} = aInferred_blind;
    all_tauRiseInferred_blind{dataset} = tauRiseInferred_blind;
    all_tauDecayInferred_blind{dataset} = tauDecayInferred_blind; 
end;

%%
figure;
all_colors = get(gca,'colororder');
close;

datasets = [1,2,3,4,5,6,7,8,9,10];
legends = {'OGB-1 (V1)', 'OGB-1 (V1)','GCaMP6s','OGB-1 (retina)',...
    'GCaMP6s', 'GCaMP5k','GCaMP6f','GCaMP6s','jRCaMP1a','jRGECO1a'};

subset_leg = [1,4,6,7,8,9,10];
colors = [4,4,3,5,3,2,1,3,6,7];
markerstyle = ['o','o','s','o','s','o','o','o','o','o'];


% datasets = [1,2,4,6,7,8,9,10];
% legends = {'OGB-1 (V1)', 'OGB-1 (V1)','OGB-1 (retina)',...
%      'GCaMP5k','GCaMP6f','GCaMP6s','jRCaMP1a','jRGECO1a'};
% 
% subset_leg = [1,3,4,5,6,7,8];
% colors = [4,4,5,2,1,3,6,7];


figPosition = [300,300,800,800];
markerSize = 50;



all_tauRise_ground_truth = cat(1,all_tauRiseInferred_groundTruth{datasets});
all_tauRise_blind= cat(1,all_tauRiseInferred_blind{datasets});

tmp = all_tauRise_blind<1.0; % Remove a weird outlier with no spikes.


all_tauRise_ground_truth = all_tauRise_ground_truth(tmp);
all_tauRise_blind = all_tauRise_blind(tmp);

correl_tauR = corrcoef(all_tauRise_ground_truth,all_tauRise_blind); correl_tauR = correl_tauR(1,2);


all_tauDecay_ground_truth = cat(1,all_tauDecayInferred_groundTruth{datasets});
all_tauDecay_blind= cat(1,all_tauDecayInferred_blind{datasets});

all_tauDecay_ground_truth = all_tauDecay_ground_truth(tmp);
all_tauDecay_blind = all_tauDecay_blind(tmp);

correl_tauD = corrcoef(all_tauDecay_ground_truth,all_tauDecay_blind); correl_tauD = correl_tauD(1,2);


all_a_ground_truth = cat(1,all_aInferred_groundTruth{datasets});
all_a_blind = cat(1,all_aInferred_blind{datasets});

all_a_ground_truth  = all_a_ground_truth(tmp);
all_a_blind  = all_a_blind(tmp);

correl_a = corrcoef(all_a_ground_truth,all_a_blind); correl_a = correl_a(1,2);



%%

figure('Position',figPosition);
hold on;
scats = cell(length(datasets),1);
for j = 1:length(datasets);
    scats{j} = scatter(all_aInferred_groundTruth{datasets(j)}, all_aInferred_blind{datasets(j)},markerSize,all_colors(colors(j),:),'Marker',markerstyle(datasets(j)),'MarkerFaceColor', all_colors(colors(j),:) );
end;

mini = 0;
maxi = 4.0;

plot([mini,maxi],[mini,maxi],'LineWidth',2.5,'Color','Black');
legend([scats{subset_leg(:)}], legends{subset_leg(:)},  'Location','SouthEast','NumColumns',2);
legend boxoff;
xlim([mini,maxi])
ylim([mini,maxi])
xlabel('Ground Truth a');
ylabel('Blind inference a');
title(sprintf('Distribution of a: correlation %.2f (0.29)',correl_a));
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('fig_spikeFinder_evaluateinference_a2.fig');
export_fig('fig_spikeFinder_evaluateinference_a2.png','-dpng');
close;

%%

figure('Position',figPosition);
hold on;
scats = cell(length(datasets),1);
for j = 1:length(datasets);
    scats{j} = scatter(all_tauRiseInferred_groundTruth{datasets(j)}, all_tauRiseInferred_blind{datasets(j)},markerSize,all_colors(colors(j),:),'Marker',markerstyle(datasets(j)),'MarkerFaceColor', all_colors(colors(j),:));
end;

mini = 0;
maxi = 0.3;

plot([mini,maxi],[mini,maxi],'LineWidth',2.5,'Color','Black');
legend([scats{subset_leg(:)}], legends{subset_leg(:)},  'Location','SouthEast','NumColumns',2);
legend boxoff;
xlim([mini,maxi])
ylim([mini,maxi])
xlabel('Ground Truth \tau_r');
ylabel('Blind inference \tau_r');
title(sprintf('Distribution of \\tau_r: correlation %.2f (0.62)',correl_tauR));
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('fig_spikeFinder_evaluateinference_taur2.fig');
export_fig('fig_spikeFinder_evaluateinference_taur2.png','-dpng');
close;

%%

figure('Position',figPosition);
hold on;
scats = cell(length(datasets),1);
for j = 1:length(datasets);
    scats{j} = scatter(all_tauDecayInferred_groundTruth{datasets(j)}, all_tauDecayInferred_blind{datasets(j)},markerSize,all_colors(colors(j),:),'Marker',markerstyle(datasets(j)),'MarkerFaceColor', all_colors(colors(j),:));
end;

mini = 0.25;
maxi = 1.4;

plot([mini,maxi],[mini,maxi],'LineWidth',2.5,'Color','Black');
legend([scats{subset_leg(:)}], legends{subset_leg(:)},  'Location','SouthEast','NumColumns',2);
legend boxoff;
xlim([mini,maxi])
ylim([mini,maxi])
xlabel('Ground Truth \tau_d');
ylabel('Blind inference \tau_d');
title(sprintf('Distribution of \\tau_d: correlation %.2f (0.69)',correl_tauD));
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('fig_spikeFinder_evaluateinference_taud2.fig');
export_fig('fig_spikeFinder_evaluateinference_taud2.png','-dpng');
close;
