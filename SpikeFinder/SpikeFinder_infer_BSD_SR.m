addpath(genpath('../'));    
load('ground_truth_time_constants.mat');

folder = 'predictions_SR/';
mkdir(folder);

q = 0.15;

tauBaselines = [60, 10, 60, 30, 40, 40, 10, 60, 30, 20]; % Values determined on train set.
superResolution = [1,1,1,1,1,2,2,2,6,3];

for k = 6:10; 
    load(sprintf('Datasets/%d.train.mat',k) );  
    T = double(T);    
    P = struct;
    P.tauRise = tauRiseGroundTruth(k); % Initial value based on ground truth inference.
    P.tauDecay = tauDecayGroundTruth(k); % Initial value based on ground truth inference.
    
    Ninf = cell(1,nNeurons);
    Cinf = cell(1,nNeurons);
    Oinf = cell(1,nNeurons);
    pAlg = cell(1,nNeurons);
    pPhys = cell(1,nNeurons);
    for i=1:nNeurons;
        O = struct;
        O.dt = dt(i);
        O.nNeurons = 1;
        O.Time = T(i);
        O.iterations = 200; 
        O.superResolution = superResolution(k);
        O.tauRiseMin = tauRiseMin(k);
        O.tauDecayMin = tauDecayMin(k);
        O.tauRiseMax = tauRiseMax(k);
        O.tauDecayMax = tauDecayMax(k);
        tmpF = normalize_remove_baseline(F{i},q,tauBaselines(k),O.dt);
        [Ninf{i},Cinf{i},pAlg{i},pPhys{i},Oinf{i}] = BSD(tmpF,O,P);
        Ninf{i} = Ninf{i}/pAlg{i}.a;
    end;
    save(sprintf('%s/%d.train.mat',folder,k ) );
end;

