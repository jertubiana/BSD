%% Script for generating figure 8.
% To execute this script, need to download the original SpikeFinder csv
% files. They contain the same spike trains, but binned at higher
% resolution (100Hz).

addpath(genpath('../'));

algorithms = {'adaptive BSD (blind)','adaptive BSD + SR'};

folders = {
    '../SpikeFinder/predictions_groundtruth/',...
    '../SpikeFinder/predictions_SR/'...    
    };

nAlg = length(algorithms);

datasets = [6,7,8,9,10];
nDatasets = length(datasets);
DeltaT = [-50:50] * 0.01;


%% Compute empirical PSFs.


all_empirical_PSFs = cell(nDatasets,nAlg);
all_mean_PSFs = cell(nDatasets,nAlg);

tauBaselines = [60, 10, 60, 30, 40, 40, 10, 60, 30, 20]; % Determined by cross-validation.
time_offsets = [0,-2,2,1,-9, 0,0,0,0,0]; % Offset between fluorescence and electrophysiological recording (in units of 0.01s). Determined by cross-validation.
dts = [1/40, 1/12 ,1/60,1/8 ,1/60,1/50, 1/60, 1/60, 1/15 , 1/30];

for k =1:nDatasets;
    dataset = datasets(k);
    Spikes = csvread( sprintf('../SpikeFinder/RawDatasets/%d.train.spikes.csv',k) );
    for u = 1 : nAlg;
        folder = folders{u};
        load(sprintf('%s/%d.train.mat',folder,dataset) );
        empirical_PSF = zeros(nNeurons,101);
        numSpikes = zeros(1,nNeurons);
        superResolution = Oinf{1}.superResolution;
        for i = 1:nNeurons;
            N = Spikes(:,i);
            % First, resample the inferred spikes at 100Hz.
            time_spikes = 0.01 * [0:size(N,1)];
            
            time_fluorescence = Time{i}; % Fluorescence measurement times.
            time_fluorescence = time_fluorescence - time_fluorescence(1);
            time_fluorescence = time_fluorescence + time_offsets(dataset) * 0.01; % For datasets 1-5, add a dataset-dependent offset.
            
            if superResolution == 1
                time_spikes_inferred = time_fluorescence;
            else
                time_spikes_inferred = interp1([1:size(N,1)+1], time_fluorescence, [1:size(N,1)*superResolution+1]/superResolution, 'linear','extrap');
            end;
            
            Ninferred = interp1(time_spikes_inferred,Ninf{i},time_spikes,'next','pchip',nan);
                        
            indices = (~isnan(N) ) & (~isnan(Ninferred) );
            
            empirical_PSF(i,:) = estimatePSF(N(indices), Ninf(indices), 50);
            empirical_PSF(i,:) = empirical_PSF(i,:)/max(empirical_PSF(i,:));
            numSpikes(i) = sum(N(indices));
        end;
        mean_PSF = (numSpikes*empirical_PSF)/sum(numSpikes);
        all_empirical_PSFs{k,u} = empirical_PSF;
        all_mean_PSFs{k,u} = mean_PSF;
    end;    
end;


%% Compute predicted PSFs.

load('../SpikeFinder/ground_truth_time_constants.mat');

PGroundTruth = cell(nDatasets,1);
OGroundTruth = cell(nDatasets,1);

all_predicted_PSFs = zeros(nDatasets,101);

for k = 1:nDatasets
    dataset = datasets(k);
    PGroundTruth{k} = struct;
    PGroundTruth{k}.SNR = median(all_aInferred{dataset}/all_sigmaInferred{dataset});
    PGroundTruth{k}.tauRise= tauRiseGroundTruth{dataset};
    PGroundTruth{k}.tauDecay= tauDecayGroundTruth{dataset};
    OGroundTruth{k} = struct;
    OGroundTruth{k}.dt = dts(dataset);
    OGroundTruth{k}.nNeurons = 1;
    OGroundTruth{k}.Nsim = 1e5;
    all_predicted_PSFs(k,:) = BSD_theoretical_accuracy(PGroundTruth,OGroundTruth);
end;



save('data_figure8.mat');
%% Plot
addpath(genpath('../'));
load('data_figure8.mat');
calciumIndicators = {'GCaMP5k','GCaMP6f','GCaMP6s','jRCAMP1a','jRGECO1a'};
BSD_functions;
s = 1;

figure;
all_default_colors = get(gca,'colororder');
close;

for k = 1:nDatasets;
    dataset = datasets(k);
    DeltaT_sampling = dts(dataset) * [-50:50];
    DeltaT_kernel = dts(dataset) * [-50:0.1:50];
    DeltaT_PSF = [-50:50] * 0.01;
    sigma = 1/PGroundTruth{k}.SNR;
    tauRise = PGroundTruth{k}.tauRise;
    tauDecay = PGroundTruth{k}.tauDecay;
    calciumIndicator = calciumIndicators{k};
    
    kernel = Kernel(tauRise,tauDecay,[-50:0.1:50],dts(dataset) );
    kernel_sampled = Kernel(tauRise,tauDecay,[-50:50],dts(dataset) );

    kernel_d = kernel - sigma;
    kernel_u = kernel + sigma;
    
    
    PSF = all_mean_PSFs{k,1};
    PSF_SR = all_mean_PSFs{k,2};
    PSF_predicted = all_predicted_PSFs(k,:);
    PSF_predicted(isnan(PSF_predicted))=0;
    
    PSF = PSF/sum(PSF);
    PSF_SR = PSF_SR/sum(PSF_SR);
    PSF_predicted = PSF_predicted /sum(PSF_predicted);    
    scale =0.25;

    X = [DeltaT_kernel,fliplr(DeltaT_kernel)];
    Y = [scale*kernel_d,fliplr(scale*kernel_u)];

    figure; 
    alpha = 0.8;
    fill(X,Y,alpha * [1,1,1]); hold on;
    h1=plot(DeltaT_PSF,PSF,'LineWidth',2); 
    h2=plot(DeltaT_PSF,PSF_SR,'LineWidth',2);
    h3=plot(DeltaT_kernel,scale * kernel);
    h4=scatter(DeltaT_sampling,scale * kernel_sampled,50,all_default_colors(4,:),'Marker','*');
    plot([0,0],[min(Y(:))-0.01,max(Y(:))],'LineWidth',2);
    h5=plot(DeltaT_PSF,PSF_predicted,'LineWidth',1,'LineStyle','--');    
    xlabel('Time (s)');
    set(gca,'Ytick',[]);
    if (k == 3) || (k == 4);
        xlim([-0.2,0.2]);
    else
        xlim([-0.1,0.1]);
    end;
    legend([h1 h2 h5 h3],'Empirical (adaptive BSD)','Empirical (adaptive BSD + SR)','Predicted','Kernel','location','best');
    title(sprintf('Point-spread function   %s',calciumIndicator) );
    legend boxoff
    set(gca,'FontSize',16);    
    
    savefig(sprintf('figure8_%s.fig',calciumIndicator));
    set(gcf,'Color','White');
    export_fig(sprintf('figure8_%s.png',calciumIndicator),'-nocrop');
    
end;


