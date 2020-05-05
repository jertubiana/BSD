%% Script for figure 2,11,12 (without constrained-oopsi).


addpath(genpath('../'));


tauRise = 0.1;
tauDecay = 0.5;

% frequencies = [1,2,5,10,20,50,100]; % values in paper.
frequencies = [1,10,20]; % quick test values.

% VecSNR = [0.1:0.1:2]; % values in paper.
VecSNR = [0.1,0.5,1.0]; % quick test values.

nu = 0.1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
% Nspikes = 5000; % vlaue in paper.
Nspikes = 100; % quick test value.


BSD_functions % load some predefined functions.

frequency_gen = 100;
dt_gen =  1/frequency_gen;
Time_gen = round(Nspikes/ (nu/frequency_gen) );

Ogen = struct; 
Ogen.Time = Time_gen;
Ogen.dt = dt_gen;
Pgen = struct;
Pgen.a = 1;
Pgen.b = 0;
Pgen.tauRise = tauRise;
Pgen.tauDecay = tauDecay;
Pgen.nu = nu;
Pgen.delta=deltaCoeff(Pgen.tauRise,Pgen.tauDecay,Ogen.dt);
Pgen.gamma=gammaCoeff(Pgen.tauRise,Pgen.tauDecay,Ogen.dt);
Pgen.eta=etaCoeff(Pgen.tauRise,Pgen.tauDecay,Ogen.dt);


N = poissrnd(Pgen.nu*Ogen.dt*ones(Ogen.Time,1)); % simulate spike train
noise = randn(Ogen.Time,1); % noise
C = Pgen.a * filter(1,[1 -Pgen.gamma-Pgen.delta,Pgen.delta],N)/Pgen.eta; % convolved spike trains.

frequencies_test = [0,1,5,10,20,100];

nFrequencies = length(frequencies);
nFrequenciesTest = length(frequencies_test);
nSNR = length(VecSNR);
nExperiments = nFrequencies * nSNR;

time_BSD = zeros(nExperiments,1);
time_oopsi = zeros(nExperiments,1);
time_nonnegative = zeros(nExperiments,1);
% time_coopsi = zeros(nExperiments,1);

correlations_BSD = zeros(nExperiments,nFrequenciesTest);
correlations_oopsi = zeros(nExperiments,nFrequenciesTest);
correlations_nonnegative = zeros(nExperiments,nFrequenciesTest);
% correlations_coopsi = zeros(nExperiments,nFrequenciesTest);



parfor experiment = 1: nExperiments
    k = 1 + floor( (experiment-1)/nSNR );
    l = 1 + mod(experiment-1, nSNR);
    
    f = frequencies(k);
    dt = 1/f;
    SNR = VecSNR(l);
    
    step = dt/dt_gen;
    F = C(step:step:end) + SNR * noise(step:step:end);
    
    O = struct;
    O.dt = dt;
    O.Time = length(F);
    O.nNeurons = 1;
    O.adaptive = 0;
    P = struct;
    P.sigma = SNR;
    P.a = 1;
    P.b = 0;
    P.nu = nu;
    P.tauRise =tauRise;
    P.tauDecay = tauDecay;
    P.delta=deltaCoeff(P.tauRise,P.tauDecay,O.dt);
    P.gamma=gammaCoeff(P.tauRise,P.tauDecay,O.dt);
    P.eta=etaCoeff(P.tauRise,P.tauDecay,O.dt);

    
    t = tic;
    Ninf_BSD =BSD( F , O , P);
    time_BSD(experiment) = toc(t);
    
    P.lambda = (P.sigma^2/(P.nu*P.a*O.dt));
    t = tic;
    Ninf_oopsi =BSD( F , O , P); % Refer to BBSD for the output.
    time_oopsi(experiment) = toc(t);
    
    P.lambda = 0;
    t = tic;
    Ninf_nonnegative =BSD( F , O , P); % Refer to BBSD for the output.
    time_nonnegative(experiment) = toc(t);    
    
%     g = [P.gamma + P.delta, -P.delta];
%     b = 0;
%     c1 = [];
%     t = tic;
%     [~,~,~,~,~,Ninf_coopsi] = constrained_foopsi(F,b,c1,g);
%     time_coopsi(experiment) = toc(t);
    
    
    for m = 1:nFrequenciesTest
        frequency_test = frequencies_test(m);
        if frequency_test ==0;
            frequency_test = f;
        else
            if frequency_test>f
                frequency_test = floor(frequency_test/f) * f;
            else
                frequency_test = f/floor(f/frequency_test) ;
            end;
        end;
        
        Nresampled = resample_spike(N,frequency_gen,frequency_test);
        
        Ninf_resampled = resample_spike(Ninf_BSD,f,frequency_test);
        correl = corrcoef(Nresampled,Ninf_resampled);
        correlations_BSD(experiment,m) = correl(1,2);
        
        Ninf_resampled = resample_spike(Ninf_oopsi,f,frequency_test);
        correl = corrcoef(Nresampled,Ninf_resampled);
        correlations_oopsi(experiment,m) = correl(1,2);
        
        
        Ninf_resampled = resample_spike(Ninf_nonnegative,f,frequency_test);   
        correl = corrcoef(Nresampled,Ninf_resampled);
        correlations_nonnegative(experiment,m) = correl(1,2);
        
        
%         Ninf_resampled = resample_spike(Ninf_coopsi,f,frequency_test);        
%         correl = corrcoef(Nresampled,Ninf_resampled);
%         correlations_coopsi(experiment,m) = correl(1,2);
    end;
end;
    
        
save('data_figure2.mat');








%%

load('data_figure2.mat');
time_BSD = reshape(time_BSD,[nFrequencies,nSNR]);
time_oopsi = reshape(time_oopsi,[nFrequencies,nSNR]);
time_nonnegative = reshape(time_nonnegative,[nFrequencies,nSNR]);
% time_coopsi = reshape(time_coopsi,[nFrequencies,nSNR]);

correlations_BSD = reshape(correlations_BSD,[nSNR,nFrequencies,nFrequenciesTest]);
correlations_oopsi = reshape(correlations_oopsi,[nSNR,nFrequencies,nFrequenciesTest]);
correlations_nonnegative = reshape(correlations_nonnegative,[nSNR,nFrequencies,nFrequenciesTest]);
% correlations_coopsi = reshape(correlations_coopsi,[nSNR,nFrequencies,nFrequenciesTest]);



for k = 1:nFrequencies
    for m = 1:nFrequenciesTest
        f = frequencies(k);
        frequency_test = frequencies_test(m);
        if frequency_test ==0;
            frequency_test = f;
        else
            if frequency_test>f
                frequency_test = floor(frequency_test/f) * f;
            else
                frequency_test = f/floor(f/frequency_test) ;
            end;
        end;

        figure('Position',[300,300,800,800]); hold on;
        p = plot(VecSNR, correlations_BSD(:,k,m),'LineWidth',3);
%         plot(VecSNR,correlations_coopsi(:,k,m),'LineWidth',3);        
        plot(VecSNR, correlations_oopsi(:,k,m),'LineWidth',3);
        plot(VecSNR,correlations_nonnegative(:,k,m),'LineWidth',3);
        if f ==1;
%             legend({'BSD','con-oopsi','oopsi','non-negative'},'Location','best','FontSize',40);
            legend({'BSD','oopsi','non-negative'},'Location','best','FontSize',40);
            
            
            legend boxoff    % Hides the legend's axes (legend border and background)
        end;

        xlim([VecSNR(1),VecSNR(end)]);
        ymin = 0; 
        ymax = ceil(correlations_BSD(1,k,m)*10)/10;
        ylim([ymin,ymax]);
        set(gca,'FontSize',35);
        xlabel('Noise to signal ratio $\frac{\sigma}{a}$','Interpreter','Latex');
        ylabel(sprintf('Correlation coefficient (%.f Hz)',frequency_test));
        title(sprintf('$f = %.f$ Hz',frequencies(k)),'Interpreter','LateX');
        set(gcf, 'Color', 'w');
    end;
end;