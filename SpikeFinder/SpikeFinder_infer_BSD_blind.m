addpath(genpath('../'));    

folder = 'predictions_blind/';
mkdir(folder);

q = 0.15;

tauBaselines = [60, 10, 60, 30, 40, 40, 10, 60, 30, 20]; % Values determined on train set.


for k = 1:10;
    load(sprintf('Datasets/%d.train.mat',k) );
    T = double(T);    
    
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
        tmpF = normalize_remove_baseline(F{i},q,tauBaselines(k),O.dt);            
        [Ninf{i},Cinf{i},pAlg{i},pPhys{i},Oinf{i}] = BSD(tmpF,O,P);
        Ninf{i} = Ninf{i}/pAlg{i}.a;
    end;
    save(sprintf('%s/%d.train.mat',folder,k) );
end;

%


for k = 1:5;
    load(sprintf('Datasets/%d.test.mat',k) );
    T = double(T);        
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
        tmpF = normalize_remove_baseline(F{i},q,tauBaselines(k),O.dt);            
        [Ninf{i},Cinf{i},pAlg{i},pPhys{i},Oinf{i}] = BSD(tmpF,O,P);
        Ninf{i} = Ninf{i}/pAlg{i}.a;
    end;
    save(sprintf('%s/%d.test.mat',folder,k) );
end;