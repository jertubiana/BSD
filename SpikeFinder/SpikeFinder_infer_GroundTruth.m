%% Ground truth parameters inference SpikeFinder.
% For each joint fluorescence/electrophysiological recording in the
% SpikeFinder dataset, infer 'ground-truth' generative model parameters
% (a,tauR,tauD,sigma) from the knowledge of the spike positions.


addpath(genpath('../'));




%% Define initial values of the kernel parameters from available literature time constants.
% See sources in main text. Because several values are excessively small, in
% practice we initialize the inference from:
% tauRise' =  tauRiseLiterature + 0.25 * tauDecayLiterature
% tauDecay' = 2 * tauDecayLiterature,

tauDecay_GCaMP5k = 0.2668;
tauDecay_GCaMP6f = 0.129;
tauDecay_GCaMP6s = 0.3804;
tauDecay_jRCaMP1a = 1.3828;
tauDecay_jRGECO1a = 0.3155;
tauDecay_OGB1 = 0.4;
tauRise_GCaMP5k = 0.0221;
tauRise_GCaMP6f = 0.0235;
tauRise_GCaMP6s = 0.0980;
tauRise_jRCaMP1a = 0.0065;
tauRise_jRGECO1a = 0.062;
tauRise_OGB1 = 0.01;


tauRiseLiterature =  ...
[tauRise_OGB1,tauRise_OGB1,tauRise_GCaMP6s,tauRise_OGB1,tauRise_GCaMP6s, ...
    tauRise_GCaMP5k,tauRise_GCaMP6f,tauRise_GCaMP6s,tauRise_jRCaMP1a,tauRise_jRGECO1a];
    
tauDecayLiterature =  ...
[tauDecay_OGB1,tauDecay_OGB1,tauDecay_GCaMP6s,tauDecay_OGB1,tauDecay_GCaMP6s, ...
    tauDecay_GCaMP5k,tauDecay_GCaMP6f,tauDecay_GCaMP6s,tauDecay_jRCaMP1a,tauDecay_jRGECO1a];


%% Infer ground truth parameters. 
% See main text for implementation details.

% Inference hyperparameters
tauBaseline = 10; % Baseline fluctuation time scale (s).
offset = 0; % Offset (in time frames), between the spikes and fluorescence.

all_tauRiseInferred = cell(10,1);
all_tauDecayInferred = cell(10,1);
all_aInferred = cell(10,1);
all_bInferred = cell(10,1);
all_sigmaInferred = cell(10,1);
tauRiseGroundTruth = zeros(10,1);
tauDecayGroundTruth = zeros(10,1);
tauRiseMin = zeros(10,1);
tauDecayMin = zeros(10,1);    
tauRiseMax = zeros(10,1);
tauDecayMax = zeros(10,1);



for k = 1:10;
    load(sprintf('../SpikeFinder/Datasets/%d.train.mat',k) );
    tauRiseInferred_ = zeros(nNeurons,1);
    tauDecayInferred_ = zeros(nNeurons,1);
    aInferred_ = zeros(nNeurons,1);
    bInferred_ = zeros(nNeurons,1);
    sigmaInferred_ = zeros(nNeurons,1);
    
    for i = 1:   nNeurons;
        Fluorescence = F{i};
        Fluorescence = (Fluorescence - median(Fluorescence) )/std(Fluorescence);
        % Normalize each recording to zero median and unit standard
        % deviation.

        TrueSpike = N{i};        
        % For datasets 3 and 5, there is a noticeable offset between the
        % fluorescence recording and the fluorescence recording. First, smooth the
        % spike train over a 10 frame window.
        
        if (k ==3) || (k==5); 
            tmp = TrueSpike;
            for t = 1+10:length(TrueSpike)-10;
                tmp(t) = max(TrueSpike(t-10:t+10));
            end;
            TrueSpike = tmp;
        end;
        
        
        tauRise0 = tauRiseLiterature(k)+0.25*tauDecayLiterature(k);
        tauDecay0 = tauDecayLiterature(k)*2;
        
        try
            [Ninf,C,aInferred_(i),bInferred_(i),sigmaInferred_(i),tauRiseInferred_(i),tauDecayInferred_(i)] = ...
                infer_generative_params_from_groundtruth(Fluorescence,TrueSpike,dt(i), tauBaseline, tauRise0, tauDecay0,offset);
        catch
            aInferred_(i) = NaN;
            tauRiseInferred_(i) = NaN;
            tauDecayInferred_(i) = NaN;
            bInferred_(i) = NaN;
            sigmaInferred_(i) = NaN;
        end
    end;
    median_a = median( aInferred_(~isnan(aInferred_) ) );
    aInferred_(isnan(aInferred_) ) = median_a;
    median_b = median( bInferred_(~isnan(bInferred_) ) );
    bInferred_(isnan(bInferred_) ) = median_b;
    median_sigma = median( sigmaInferred_(~isnan(sigmaInferred_) ) );
    sigmaInferred_(isnan(sigmaInferred_) ) = median_sigma;
    median_tauR = median( tauRiseInferred_(~isnan(tauRiseInferred_) ) );
    tauRiseInferred_(isnan(tauRiseInferred_) ) = median_tauR;
    median_tauD = median( tauDecayInferred_(~isnan(tauDecayInferred_) ) );
    tauDecayInferred_(isnan(tauDecayInferred_) ) = median_tauD;
    
    
    
    all_tauRiseInferred{k} = tauRiseInferred_;
    all_tauDecayInferred{k} = tauDecayInferred_;
    all_aInferred{k} = aInferred_;
    all_bInferred{k} = bInferred_;
    all_sigmaInferred{k} = sigmaInferred_;
    tauRiseGroundTruth(k) = median(all_tauRiseInferred{k});
    tauDecayGroundTruth(k) = median(all_tauDecayInferred{k});
    tauRiseMin(k) = min(all_tauRiseInferred{k});
    tauDecayMin(k) = min(all_tauDecayInferred{k});
    tauRiseMax(k) = max(all_tauRiseInferred{k});
    tauDecayMax(k) = max(all_tauDecayInferred{k});
    save('ground_truth_time_constants.mat','tauRiseGroundTruth','tauDecayGroundTruth','tauRiseMin','tauRiseMax','tauDecayMin','tauDecayMax','all_tauRiseInferred','all_tauDecayInferred','all_aInferred','all_bInferred','all_sigmaInferred');
end;    

save('ground_truth_time_constants.mat','tauRiseGroundTruth','tauDecayGroundTruth','tauRiseMin','tauRiseMax','tauDecayMin','tauDecayMax','all_tauRiseInferred','all_tauDecayInferred','all_aInferred','all_bInferred','all_sigmaInferred');

