function [Nconvolved , all_PSF,DeltaT]= BSD_convolve_with_PSF(N,O,P,varargin)
% Computes the predicted point-spread function, and convolve the spike
% train with the output.
% Inputs: 
% N: An array or a cell of (inferred) spike trains.
% O: A structure or cell of structures of mandatory inputs and optional algorithm options, same
% format as for main BSD file. Must contain nNeurons,Time,dt.
% P: A structure or cell of structures of inferred/known generative
% parameters. Must contain tauRise, tauDecay, (a, sigma) or SNR.
% Optional additional argument 'median' or 'single': either compute an individual PSF for each neuron, or
% compute a single PSF for all the dataset, using the median generative
% parameters. Default = median.
% Output:
% Typical usages:
% [Ninf,Cinf,Pinf,pPhysinf,Oinf]=BSD( Fluorescence , O)
% Ninfconvolved = BSD_convolve_with_PSF(Ninf,Oinf,Pinf);
% ------
% [Ninf,Cinf,Pinf,pPhysinf,Oinf]=BSD( Fluorescence , O)
% Ninfconvolved = BSD_convolve_with_PSF(Ninf,Oinf,Pinf,'single'); % Compute
% the PSF individually for each neuron.


N_is_cell = strcmp(class(N) , 'cell');
O_is_cell = strcmp(class(O) , 'cell');
P_is_cell = strcmp(class(P) ,'cell');

if length(varargin)>0
    mode = varargin{1};
else
    mode = 'median';
end;

case1 = (~N_is_cell) & (~O_is_cell) & (~P_is_cell); % Similar to the output of BSD-like.
case2 = (N_is_cell) & (O_is_cell) & (P_is_cell); % When each recording has have different length.
case3 = (N_is_cell) & (O_is_cell) & (~P_is_cell); % When each recording is different, but use same parameters for all. 
case4 = (N_is_cell) & (~O_is_cell) & (~P_is_cell); % When each recording is different, but use same parameters for all. 

display([case1,case2,case3,case4]);

if case1;
    nNeurons = O.nNeurons;
    if size(N,1) ~= nNeurons;
        need_transpose = true;
        N = N';
    else
        need_transpose = false;
    end;
        
    Nconvolved = zeros(O.nNeurons,O.Time);
elseif (case2 || case3 || case4);
    nNeurons = length(N);
    Nconvolved = cell(nNeurons,1);
end;


if strcmp(mode , 'median')
    all_dt = zeros(nNeurons,1);
    all_tauRise = zeros(nNeurons,1);
    all_tauDecay = zeros(nNeurons,1);
    all_SNR = zeros(nNeurons,1);
    
    for n = 1:nNeurons;
        if case1 || case4
            all_dt(n) = O.dt;
        elseif case2 || case3
            all_dt(n) = O{n}.dt;
        end;
        if case1 || case3 || case4
            try
                all_tauRise(n) = P.tauRise(n);
            catch
                all_tauRise(n) = P.tauRise;
            end;
            
            try
                all_tauDecay(n) = P.tauDecay(n);
            catch
                all_tauDecay(n) = P.tauDecay;
            end;
            
            if isfield(P,'SNR');
                try
                    all_SNR(n) = P.SNR(n);
                catch
                    all_SNR(n) = P.SNR;
                end;                
            else
                try
                    all_SNR(n) = P.a(n)/P.sigma(n);
                catch
                    all_SNR(n) = P.a/P.sigma;
                end;
            end;
        else
            all_tauRise(n) = P{n}.tauRise;
            all_tauDecay(n) = P{n}.tauDecay;
            if isfield(P,'SNR');
                all_SNR(n) = P{n}.SNR;
            else
                all_SNR(n) = P{n}.a/P{n}.sigma;
            end
        end;
    end;
        
    
    dt_median = median(all_dt);
    tauRise_median = median(all_tauRise);
    tauDecay_median = median(all_tauDecay);
    SNR_median = median(all_SNR);
    
    Oprime = struct;
    Oprime.dt = dt_median;
    if case1 || case4;
        if isfield(O,'z1'); Oprime.z1 = O.z1; end;
        if isfield(O,'z2'); Oprime.z2 = O.z2; end;
        if isfield(O,'delta_max'); Oprime.delta_max = O.delta_max; end;
        if isfield(O,'Nsim'); Oprime.Nsim = O.Nsim; end;
        if isfield(O,'superResolution'); Oprime.superResolution = O.superResolution; end;
        if isfield(O,'discretization'); Oprime.discretization = O.discretization; end;            
    elseif case2 || case3;
        n = 1;
        Oprime.dt = O{n}.dt;
        if isfield(O{n},'z1'); Oprime.z1 = O{n}.z1; end;
        if isfield(O{n},'z2'); Oprime.z2 = O{n}.z2; end;
        if isfield(O{n},'delta_max'); Oprime.delta_max = O{n}.delta_max; end;
        if isfield(O{n},'Nsim'); Oprime.Nsim = O{n}.Nsim; end;
        if isfield(O{n},'superResolution'); Oprime.superResolution = O{n}.superResolution; end;
        if isfield(O{n},'discretization'); Oprime.discretization = O{n}.discretization; end;
    end;
        
    Pprime = struct;
    Pprime.tauRise = tauRise_median;
    Pprime.tauDecay = tauDecay_median;
    Pprime.sigma = 1/SNR_median;
    Pprime.a = 1;
    [all_PSF,DeltaT,~,~,~,~] = BSD_theoretical_accuracy(Pprime,Oprime);
    [maxVal,location] = max(all_PSF);
    DeltaT_target = Oprime.dt * [-30:30];
    all_PSF_discrete = interp1(DeltaT-DeltaT(location),all_PSF,DeltaT_target);
    all_PSF_discrete(isnan(all_PSF_discrete))=0;

    
    
    
elseif strcmp(mode , 'single')
    all_PSF = cell(nNeurons,1);
    all_PSF_discrete = cell(nNeurons,1);
    for n =1:nNeurons;
        Oprime = struct;
        if case1 || case4;
            Oprime.dt = O.dt;
            if isfield(O,'z1'); Oprime.z1 = O.z1; end;
            if isfield(O,'z2'); Oprime.z2 = O.z2; end;
            if isfield(O,'delta_max'); Oprime.delta_max = O.delta_max; end;
            if isfield(O,'Nsim'); Oprime.Nsim = O.Nsim; end;
            if isfield(O,'superResolution'); Oprime.superResolution = O.superResolution; end;
            if isfield(O,'discretization'); Oprime.discretization = O.discretization; end;            
        elseif case2 || case3;
            Oprime.dt = O{n}.dt;
            if isfield(O{n},'z1'); Oprime.z1 = O{n}.z1; end;
            if isfield(O{n},'z2'); Oprime.z2 = O{n}.z2; end;
            if isfield(O{n},'delta_max'); Oprime.delta_max = O{n}.delta_max; end;
            if isfield(O{n},'Nsim'); Oprime.Nsim = O{n}.Nsim; end;
            if isfield(O{n},'superResolution'); Oprime.superResolution = O{n}.superResolution; end;
            if isfield(O{n},'discretization'); Oprime.discretization = O{n}.discretization; end;
        end;

        Pprime = struct;
        if case1 || case3 || case4;
            Pprime.tauRise = P.tauRise(n);
            Pprime.tauDecay= P.tauDecay(n);
            if isfield(P,'sigma') && isfield(P,'a')
                Pprime.sigma = P.sigma(n);
                Pprime.a = P.a(n);
            else
                Pprime.sigma = 1/P.SNR(n);
                Pprime.a = 1;
            end;
        elseif case2
            Pprime.tauRise = P{n}.tauRise;
            Pprime.tauDecay= P{n}.tauDecay;
            if isfield(P{n},'sigma') && isfield(P{n},'a')
                Pprime.sigma = P{n}.sigma;
                Pprime.a = P{n}.a;
            else
                Pprime.sigma = 1/P{n}.SNR;
                Pprime.a = 1;
            end;            
        end;
        [all_PSF{n},DeltaT,~,~,~,~] = BSD_theoretical_accuracy(Pprime,Oprime);
        display(all_PSF{n});
        DeltaT_target = Oprime.dt * [-30:30];
        [maxVal,location] = max(all_PSF{n}); 
        PSF_discrete = interp1(DeltaT-DeltaT(location),all_PSF{n},DeltaT_target);
        PSF_discrete(isnan(PSF_discrete))=0;
        all_PSF_discrete{n} = PSF_discrete;
    end;
end;


for n = 1:nNeurons;
    if strcmp(mode , 'single')
        PSF_discrete = all_PSF_discrete{n};
    else
        PSF_discrete = all_PSF_discrete;
    end

    if case1;
        Nconvolved(n,:) = conv(N(n,:), PSF_discrete/max(PSF_discrete),'same');
    else
        Nconvolved{n} = conv(N{n}, PSF_discrete/max(PSF_discrete),'same');
    end;
end;

if case1 && need_transpose;
    Nconvolved = Nconvolved';
end;

end