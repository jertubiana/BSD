%% Script for figures 3 and 13.
% 1) Generate a "synthetic version" of the SpikeFinder dataset, i.e. for
% each trace of the SpikeFinder set, infer its generative parameters using
% the spikes and the fluorescence, and then generate a few synthetic traces given the true spikes and generative parameters.

% 2) Infer spikes and parameters using the fluorescence traces only.
% 3) Compare inferred parameter values and ground truth values.
% 4) Evaluate sources of errors as function of generative parameters and spike distribution.


%% Generate a synthetic version of the SpikeFinder dataset.

clear all
addpath(genpath('../'));
ground_truth_params = '../SpikeFinder/ground_truth_time_constants.mat';
load(ground_truth_params);


datasets = [1,2,4,6,7,8,9,10];  % Values in paper.
% datasets = [6]; % For test.

% SNR_multiplicator = [0.25, 0.5, 1.0]; % Values in paper.
SNR_multiplicator = [1.0]; % For test



nDatasets = length(datasets);
nRepeats = 5;
nSNR = length(SNR_multiplicator);


all_Spikes = cell(nDatasets,1);
all_Fluorescence = cell(nDatasets,nRepeats,nSNR);
all_tauR_true = all_tauRiseInferred;
all_tauD_true = all_tauDecayInferred;
all_SNR = cell(nDatasets,1);
all_numSpikes = cell(nDatasets,1);
all_FiringRate = cell(nDatasets,1);
all_Autocorrelation_Width = cell(nDatasets,1);
all_dt = cell(nDatasets,1);
all_nNeurons = zeros(nDatasets,1);


for i = 1:nDatasets
    dataset = datasets(i);
    load(sprintf('../SpikeFinder/Datasets/%d.train.mat',dataset) );
    all_nNeurons(i) = nNeurons;    
    Spikes = cell(nNeurons,1);
    numSpikes = zeros(nNeurons,1);
    FiringRate = zeros(nNeurons,1);
    Autocorrelation_Width = zeros(nNeurons,2);

    for j = 1:nNeurons;
        Spikes{j} = double(N{j});
        numSpikes(j) = sum(Spikes{j});
        FiringRate(j) = numSpikes(j) / (dt(j) * length(Spikes{j}));
        if numSpikes(j)>0;
            Autocorrelation_Width(j,1) = get_autocorrelation_width( Spikes{j} );
            Autocorrelation_Width(j,2) = ( Autocorrelation_Width(j,1) * dt(j) )/(all_tauD_true{dataset}(j) );
        end;
    end;
    

    all_Spikes{i} = Spikes;
    all_SNR{i} = all_aInferred{dataset}./all_sigmaInferred{dataset};
    all_numSpikes{i} = numSpikes;
    all_FiringRate{i} = FiringRate;
    all_Autocorrelation_Width{i} = Autocorrelation_Width;
    all_dt{i} = dt;


    for k=1:nRepeats;
        for l = 1:nSNR;
            Fluorescence = cell(nNeurons,1);
            for j = 1:nNeurons;
                P = struct;
                P.tauRise = all_tauR_true{dataset}(j);
                P.tauDecay = all_tauD_true{dataset}(j);
                P.a = 1;
                P.b = 0;
                P.sigma = 1/(all_SNR{i}(j) * SNR_multiplicator(l) ) ;
                O = struct;
                O.Time = length(all_Spikes{i}{j});
                O.nNeurons = 1;
                O.dt = all_dt{i}(j);
                [~,~,Fluorescence{j}] = BSD_generate_synthetic_signal(P,O, all_Spikes{i}{j} );
            end;
            all_Fluorescence{i,k,l} = Fluorescence;
        end;
    end;
end;


save('synthetic_dataset_from_SpikeFinder.mat', 'all_Spikes', 'all_Fluorescence', 'all_tauR_true', 'all_tauD_true', ...
'all_SNR', 'all_numSpikes', 'all_FiringRate', 'all_Autocorrelation_Width', ...
'all_dt', 'datasets','all_nNeurons','SNR_multiplicator','nDatasets','nRepeats','nSNR');

%% Infer spikes and parameters using the fluorescence traces only.

clear all
addpath(genpath('../'));
load('synthetic_dataset_from_SpikeFinder.mat');


all_tauR_inferred = cell(nDatasets,nRepeats,nSNR);
all_tauD_inferred = cell(nDatasets,nRepeats,nSNR);

for i = 1:nDatasets;
    nNeurons = all_nNeurons(i);
    for k = 1:nRepeats;
        for l = 1:nSNR;
            tauR_inferred = zeros(nNeurons,2);
            tauD_inferred = zeros(nNeurons,2);
            for j = 1:nNeurons;
                F = all_Fluorescence{i,k,l}{j};
                O = struct;
                O.nNeurons = 1;
                O.Time = length(F);
                O.dt = all_dt{i}(j);
                O.tauRiseMax = 0.5;
                O.tauDecayMax = 3;
                
                O.iterations = 200; % Infer with iterative refinement.
                O.adaptive = 1;
                [~,~,pAlg,~,~] = BSD(F,O);
                tauR_inferred(j,1) = pAlg.tauRise;
                tauD_inferred(j,1) = pAlg.tauDecay;
                
                
                O.iterations = 0; % Infer with only one iteration.
                O.adaptive = 0;
                [~,~,pAlg,~,~] = BSD(F,O);
                tauR_inferred(j,2) = pAlg.tauRise;
                tauD_inferred(j,2) = pAlg.tauDecay; 
                
            end;
            all_tauR_inferred{i,k,l} = tauR_inferred;
            all_tauD_inferred{i,k,l} = tauD_inferred;
        end;
    end;
end;


save('results_synthetic_dataset_from_SpikeFinder.mat', 'all_tauR_true', 'all_tauD_true', ...
'all_tauR_inferred', 'all_tauD_inferred', 'all_SNR', 'all_numSpikes', 'all_FiringRate', 'all_Autocorrelation_Width', ...
'all_dt', 'datasets','SNR_multiplicator','nDatasets','nRepeats','nSNR','all_nNeurons');


%% Compare inferred parameter values and ground truth values.
addpath(genpath('../../BSD'))
BSD_functions;
load('results_synthetic_dataset_from_SpikeFinder.mat');

% Define colors for each indicator.
figure;
all_default_colors = get(gca,'colororder');
close;

all_datasets = [1,2,3,4,5,6,7,8,9,10];
all_legends = {'OGB-1 (V1)', 'OGB-1 (V1)','GCaMP6s','OGB-1 (retina)',...
    'GCaMP6s', 'GCaMP5k','GCaMP6f','GCaMP6s','jRCaMP1a','jRGECO1a'};
all_colors = [4,4,3,5,3,2,1,3,6,7];


% Regroup all experiments into a single array of results.
all_tauRise_ground_truth_flat = zeros(sum(all_nNeurons) * nSNR * nRepeats, 1);
all_tauRise_inferred_flat = zeros(sum(all_nNeurons) * nSNR * nRepeats, 2);
all_tauDecay_ground_truth_flat = zeros(sum(all_nNeurons) * nSNR * nRepeats, 1);
all_tauDecay_inferred_flat = zeros(sum(all_nNeurons) * nSNR * nRepeats, 2);
all_SNR_flat = zeros(sum(all_nNeurons) * nSNR * nRepeats, 1);
all_dt_flat = zeros(sum(all_nNeurons) * nSNR * nRepeats, 1);
all_numSpikes_flat = zeros(sum(all_nNeurons) * nSNR * nRepeats, 1);
all_FiringRate_flat = zeros(sum(all_nNeurons) * nSNR * nRepeats, 1);
all_Autocorrelation_Width_flat = zeros(sum(all_nNeurons) * nSNR * nRepeats, 2);
all_Dataset_flat =  zeros(sum(all_nNeurons) * nSNR * nRepeats, 1);


index = 1;
for i = 1:nDatasets;
    dataset = datasets(i);
    nNeurons = all_nNeurons(i);
    for j = 1:nRepeats;
        for k = 1:nSNR;
            for l = 1:nNeurons;
                if all_numSpikes{i}(l)>0; % Only keep experiments with spikes.
                    all_tauRise_ground_truth_flat(index) = all_tauR_true{dataset}(l);
                    all_tauRise_inferred_flat(index,:) = all_tauR_inferred{i,j,k}(l,:);
                    all_tauDecay_ground_truth_flat(index) = all_tauD_true{dataset}(l);
                    all_tauDecay_inferred_flat(index,:) = all_tauD_inferred{i,j,k}(l,:);
                    all_SNR_flat(index) = all_SNR{i}(l) * SNR_multiplicator(k);
                    all_dt_flat(index) = all_dt{i}(l);
                    all_numSpikes_flat(index) = all_numSpikes{i}(l);
                    all_FiringRate_flat(index) = all_FiringRate{i}(l);
                    all_Autocorrelation_Width_flat(index,:) = all_Autocorrelation_Width{i}(l,:);
                    all_Dataset_flat(index) = dataset;
                    index = index+1;
                end;
            end;
        end;
    end;
end;

% Remove the zero entries.
all_tauRise_ground_truth_flat = all_tauRise_ground_truth_flat(1:index-1);
all_tauRise_inferred_flat = all_tauRise_inferred_flat(1:index-1,:);
all_tauDecay_ground_truth_flat = all_tauDecay_ground_truth_flat(1:index-1);
all_tauDecay_inferred_flat = all_tauDecay_inferred_flat(1:index-1,:);
all_SNR_flat = all_SNR_flat(1:index-1);
all_dt_flat = all_dt_flat(1:index-1);
all_numSpikes_flat = all_numSpikes_flat(1:index-1);
all_FiringRate_flat = all_FiringRate_flat(1:index-1);
all_Autocorrelation_Width_flat = all_Autocorrelation_Width_flat(1:index-1,:);
all_Dataset_flat = all_Dataset_flat(1:index-1);


% Compute the effective SNR sigma/a \|K\|.
all_SNR_eff_flat = zeros(size(all_SNR_flat));
for u = 1:length(all_SNR_flat);
    all_SNR_eff_flat(u) = all_SNR_flat(u)*normKernelCoeff(all_tauRise_ground_truth_flat(u),all_tauDecay_ground_truth_flat(u),all_dt_flat(u));
end;


% Compute correlation coefficients with and without iterative refinements.
correl_tauR = corrcoef(all_tauRise_ground_truth_flat,all_tauRise_inferred_flat(:,1)); correl_tauR = correl_tauR(1,2);
correl_tauD = corrcoef(all_tauDecay_ground_truth_flat,all_tauDecay_inferred_flat(:,1)); correl_tauD = correl_tauD(1,2);
correl_tauR2 = corrcoef(all_tauRise_ground_truth_flat,all_tauRise_inferred_flat(:,2)); correl_tauR2 = correl_tauR2(1,2);
correl_tauD2 = corrcoef(all_tauDecay_ground_truth_flat,all_tauDecay_inferred_flat(:,2)); correl_tauD2 = correl_tauD2(1,2);


% Compute the overlap between the ground truth and the inferred kernel
% (Higher is better).
overlap = @(tauRise1,tauDecay1,tauRise2,tauDecay2, dt) ...
(...
1.0./(1-exp(-dt./tauDecay1-dt./tauDecay2) ) + 1./(1-exp(-dt./tauRise1-dt./tauRise2) ) ...
-1./(1-exp(-dt./tauDecay1-dt./tauRise2) ) - 1./(1-exp(-dt./tauRise1-dt./tauDecay2) ) );

normedoverlap = @(tauRise1,tauDecay1,tauRise2,tauDecay2, dt) overlap(tauRise1,tauDecay1,tauRise2,tauDecay2,dt)./sqrt( overlap(tauRise1,tauDecay1,tauRise1,tauDecay1,dt) .* overlap(tauRise2,tauDecay2,tauRise2,tauDecay2,dt)  );

all_overlaps_flat = normedoverlap(all_tauRise_ground_truth_flat,all_tauDecay_ground_truth_flat,...
              all_tauRise_inferred_flat(:,1), all_tauDecay_inferred_flat(:,1),...
              all_dt_flat);


          
all_overlaps_flat = sqrt(1-all_overlaps_flat.^2);    



%% Figure 3 upper right panel.
figPosition = [300,300,800,800];
markerSize = 50;


figure('Position',figPosition);
hold on;
scats = cell(nDatasets,1);
for i = 1:nDatasets;
    dataset = datasets(i);
    subset1 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat>20);
    subset2 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat<=20);  
    
    scats{i} = scatter(all_tauRise_ground_truth_flat(subset1), all_tauRise_inferred_flat(subset1,1)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','o','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );

    scatter(all_tauRise_ground_truth_flat(subset2), all_tauRise_inferred_flat(subset2,1)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','x','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );


end;

mini = 0;
maxi = 0.3;

plot([mini,maxi],[mini,maxi],'LineWidth',2.5,'Color','Black');
legend([scats{2:end}], all_legends{datasets(2:end)},  'Location','SouthEast','NumColumns',2);
legend boxoff
xlim([mini,maxi])
% ylim([mini,maxi])
xlabel('Ground Truth \tau_r');
ylabel('Blind inference \tau_r (iterative)');
title(sprintf('Distribution of \\tau_r: correlation %.2f',correl_tauR));
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('figure_3b.fig');
export_fig('figure_3b.png','-dpng');
% close;

%% Figure 3 lower left panel.

figPosition = [300,300,800,800];
markerSize = 50;


figure('Position',figPosition);
hold on;
scats = cell(nDatasets,1);
for i = 1:nDatasets;
    dataset = datasets(i);
    subset1 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat>20);
    subset2 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat<=20);  
        
    scats{i} = scatter(all_tauRise_ground_truth_flat(subset1), all_tauRise_inferred_flat(subset1,2)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','o','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );

    scatter(all_tauRise_ground_truth_flat(subset2), all_tauRise_inferred_flat(subset2,2)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','x','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );



end;

mini = 0;
maxi = 0.3;

plot([mini,maxi],[mini,maxi],'LineWidth',2.5,'Color','Black');
legend([scats{2:end}], all_legends{datasets(2:end)},  'Location','SouthEast','NumColumns',2);
legend boxoff
xlim([mini,maxi])
% ylim([mini,maxi])
xlabel('Ground Truth \tau_r');
ylabel('Blind inference \tau_r (initial)');
title(sprintf('Distribution of \\tau_r: correlation %.2f',correl_tauR2));
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('figure_3c.fig');
export_fig('figure_3c.png','-dpng');
% close;


%% Figure 3 lower right panel.



figPosition = [300,300,800,800];
markerSize = 50;


figure('Position',figPosition);
hold on;
scats = cell(nDatasets,1);
for i = 1:nDatasets;
    dataset = datasets(i);
    subset1 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat>20);
    subset2 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat<=20);        
    
    
    scats{i} = scatter(all_tauDecay_ground_truth_flat(subset1), all_tauDecay_inferred_flat(subset1,1)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','o','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );

    scatter(all_tauDecay_ground_truth_flat(subset2), all_tauDecay_inferred_flat(subset2,1)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','x','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );


end;

mini = 0.25;
maxi = 1.4;

plot([mini,maxi],[mini,maxi],'LineWidth',2.5,'Color','Black');
legend([scats{2:end}], all_legends{datasets(2:end)},  'Location','SouthEast','NumColumns',2);
legend boxoff
xlim([mini,maxi])
% ylim([mini,maxi])
xlabel('Ground Truth \tau_d');
ylabel('Blind inference \tau_d (iterative)');
title(sprintf('Distribution of \\tau_d: correlation %.2f',correl_tauD));
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('figure_3d.fig');
export_fig('figure_3d.png','-dpng');
close;

%% Figure 3 upper left panel.

figPosition = [300,300,800,800];
markerSize = 50;



figure('Position',figPosition);
hold on;
scats = cell(nDatasets,1);
for i = 1:nDatasets;
    dataset = datasets(i);
    subset1 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat>20);
    subset2 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat<=20);  

    scats{i} = scatter(all_tauDecay_ground_truth_flat(subset1), all_tauDecay_inferred_flat(subset1,2)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','o','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );


    scatter(all_tauDecay_ground_truth_flat(subset2), all_tauDecay_inferred_flat(subset2,2)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','x','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );


end;

mini = 0.25;
maxi = 1.4;

plot([mini,maxi],[mini,maxi],'LineWidth',2.5,'Color','Black');
legend([scats{2:end}], all_legends{datasets(2:end)},  'Location','SouthEast','NumColumns',2);
legend boxoff
xlim([mini,maxi])
% ylim([mini,maxi])
xlabel('Ground Truth \tau_d');
ylabel('Blind inference \tau_d (initial)');
title(sprintf('Distribution of \\tau_d: correlation %.2f',correl_tauD2));
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('figure_3a.fig');
export_fig('figure_3a.png','-dpng');
close;


%% Overlap vs SNR, not shown in article.

figPosition = [300,300,800,800];
markerSize = 50;


figure('Position',figPosition);
hold on;
scats = cell(nDatasets,1);
for i = 1:nDatasets;
    dataset = datasets(i);
    subset1 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat>20);
    subset2 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat<=20);  
    
    scats{i} = scatter(all_SNR_flat(subset1), all_overlaps_flat(subset1)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','o','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );

    scatter(all_SNR_flat(subset2), all_overlaps_flat(subset2)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','x','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );


end;

legend([scats{2:end}], all_legends{datasets(2:end)},  'Location','best','NumColumns',2);
legend boxoff
xlabel('Signal-to-noise ratio $\frac{a}{\sigma}$','Interpreter','Latex');
ylabel('Kernel mismatch (normalized)');
title('Kernel mismatch vs SNR');
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('figure_spikefinder2synthetic_kernel_inference_mismatch_SNR.fig');
export_fig('figure_spikefinder2synthetic_kernel_inference_mismatch_SNR.png','-dpng');
% close;


%% Overlap vs effective SNR, Fig13b


figPosition = [300,300,800,800];
markerSize = 50;

% 1: nSpikes < 20.
% 2: 

figure('Position',figPosition);
hold on;
scats = cell(nDatasets,1);
for i = 1:nDatasets;
    dataset = datasets(i);
    subset1 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat>20);
    subset2 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat<=20);    
    
    scats{i} = scatter(all_SNR_eff_flat(subset1), all_overlaps_flat(subset1)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','o','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );

    scatter(all_SNR_eff_flat(subset2), all_overlaps_flat(subset2)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','x','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );


end;

legend([scats{2:end}], all_legends{datasets(2:end)},  'Location','best','NumColumns',2);
legend boxoff
xlabel('Effective Signal-to-noise ratio $\frac{a \| K\|}{\sigma}$','Interpreter','Latex');
ylabel('Kernel mismatch (normalized)');
% title('Kernel mismatch vs effective SNR');
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('figure_spikefinder2synthetic_kernel_inference_mismatch_SNReff.fig');
export_fig('figure_spikefinder2synthetic_kernel_inference_mismatch_SNReff.png','-dpng');
% close;


%% Overlap vs firing rate, not shown.

figPosition = [300,300,800,800];
markerSize = 50;


figure('Position',figPosition);
hold on;
scats = cell(nDatasets,1);
for i = 1:nDatasets;
    dataset = datasets(i);
    
    subset1 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat>20);
    subset2 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat<=20);    
    
    scats{i} = scatter(all_FiringRate_flat(subset1), all_overlaps_flat(subset1)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','o','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );

    scatter(all_FiringRate_flat(subset2), all_overlaps_flat(subset2)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','x','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );


end;

legend([scats{2:end}], all_legends{datasets(2:end)},  'Location','SouthEast','NumColumns',2);
legend boxoff
xlabel('Firing Rate (Hz)');
ylabel('Kernel mismatch (normalized)');
title('Kernel mismatch vs firing rate');
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('figure_spikefinder2synthetic_kernel_inference_mismatch_firingrate.fig');
export_fig('figure_spikefinder2synthetic_kernel_inference_mismatch_firingrate.png','-dpng');
% close;




%% Overlap vs number of spikes, not shown.

figPosition = [300,300,800,800];
markerSize = 50;


figure('Position',figPosition);
hold on;
scats = cell(nDatasets,1);
for i = 1:nDatasets;
    dataset = datasets(i);
    subset = (all_Dataset_flat == dataset );
    scats{i} = scatter(all_numSpikes_flat(subset), all_overlaps_flat(subset)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','o','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );
end;

legend([scats{2:end}], all_legends{datasets(2:end)},  'Location','best','NumColumns',2);
legend boxoff
xlabel('Number of spikes');
ylabel('Kernel mismatch (normalized)');
title('Kernel mismatch vs number of spikes');
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('figure_spikefinder2synthetic_kernel_inference_mismatch_numspikes.fig');
export_fig('figure_spikefinder2synthetic_kernel_inference_mismatch_numspikes.png','-dpng');
% close;

%% Overlap vs spike autocorrelation width, Fig 13a.

figPosition = [300,300,800,800];
markerSize = 50;


figure('Position',figPosition);
hold on;
scats = cell(nDatasets,1);
for i = 1:nDatasets;
    dataset = datasets(i);

    subset1 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat>20);
    subset2 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat<=20);        
    
    scats{i} = scatter(all_Autocorrelation_Width_flat(subset1,1), all_overlaps_flat(subset1)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','o','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );

    scatter(all_Autocorrelation_Width_flat(subset2,1), all_overlaps_flat(subset2)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','x','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );


end;

legend([scats{2:end}], all_legends{datasets(2:end)},  'Location','best','NumColumns',2);
legend boxoff
xlabel('Spike autocorrelation width (time bins)');
ylabel('Kernel mismatch (normalized)');
title('Kernel mismatch vs Spike autocorrelation width');
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('figure_spikefinder2synthetic_kernel_inference_mismatch_autocorrelation.fig');
export_fig('figure_spikefinder2synthetic_kernel_inference_mismatch_autocorrelation.png','-dpng');
% close;


%% Overlap vs spike autocorrelation width for initial inference, not shown.

figPosition = [300,300,800,800];
markerSize = 50;


maxival = 1;

figure('Position',figPosition);
hold on;
scats = cell(nDatasets,1);
for i = 1:nDatasets;
    dataset = datasets(i);

    subset1 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat>20);
    subset2 = (all_Dataset_flat == dataset ) & (all_numSpikes_flat<=20);        
    
    
    scats{i} = scatter(min(all_Autocorrelation_Width_flat(subset1,2),maxival), all_overlaps_flat(subset1)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','o','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );

    scatter(min(all_Autocorrelation_Width_flat(subset2,2),maxival), all_overlaps_flat(subset2)...
    ,markerSize,all_default_colors(all_colors(dataset),:),'Marker','x','MarkerFaceColor',  all_default_colors(all_colors(dataset),:) );

end;

legend([scats{2:end}], all_legends{datasets(2:end)},  'Location','best','NumColumns',2);
legend boxoff
xlabel('Spike autocorrelation width (normed)');
xlim([0,maxival]);
ylabel('Kernel mismatch (normalized)');
% title('Kernel mismatch vs Spike autocorrelation width (normed)');
set(gca,'FontSize',25);
set(gcf,'Color','White');
savefig('figure_spikefinder2synthetic_kernel_inference_mismatch_autocorrelation2.fig');
export_fig('figure_spikefinder2synthetic_kernel_inference_mismatch_autocorrelation2.png','-dpng');
% close;


