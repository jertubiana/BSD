%% BSD super-resolution.
% An introductory script to test BSD with super-resolution on synthetic fluorescence data.

%% Generate a fluorescence signal
clear all;
addpath(genpath('../'));
O = struct; % Struct of experimental conditions & decoding options.
O.Time = 1e4; % Number of time frames.
O.nNeurons =1; % Number of neurons.
f = 100; % Frequency of acquisition. (Hz)
O.dt = 1/f; % interval duration. (s)

P = struct; % Struct of generative model properties.
P.tauRise = 0.1; % Fluorescence raise time (s)
P.tauDecay = 0.5; % Fluorescence decay time (s)
P.nu = 0.5; % Neuron firing rate (Hz).
P.a = 1; % Spike amplitude
P.b = 0; % Baseline position
P.sigma = 0.2; % Noise standard deviation.

T= [1:O.Time]/f; % Convention: measurement times are i * \Delta t, with i>=1.
[N,C,F] = BSD_generate_synthetic_signal(P,O);
SR = 5;
downF = F(SR:SR:end); 
downT = T(SR:SR:end); % After down-sampling, the measurement times are the j * SR * \Delta t, with j>=1.


figure; plot(T,F); hold on; plot(downT,downF); plot(T,N);
xlabel('Time (s)');
legend('Fluorescence','Low resolution Fluorescence','Spikes');
set(gca,'FontSize',16)


%% Perform SR inference with known generative parameters.
Oalg = struct; % Struct of experimental conditions & decoding options.
Oalg.Time = length(downF); % Number of time frames.
Oalg.dt = O.dt*SR; % interval duration. (s)
Oalg.nNeurons = O.nNeurons; % Number of neurons.
Oalg.adaptive = 0; % Not adaptive. Will use provided values for parameters, and estimate the unknown ones.
Oalg.superResolution = SR; % Attempt to reconstruct at the original signal frequency.

Palg=P; % Use known generative parameters

tic;
[Ninf,Cinf,~,Pphys]=BSD( downF , Oalg , Palg);
toc;
binNinf = Ninf > Pphys.threshold;


figure(1); 
h1 = subplot(2,1,1);
plot(downT,downF,'Color',[0 0 0] + 0.5); hold on; 
plot(T,Cinf,'green'); legend('DF/F', 'Denoised DF/F');
set(gca,'FontSize',16);
h2 = subplot(2,1,2);
plot(T,N,'black'); hold on;
plot(T,Ninf,'green'); plot(T,binNinf*0.5); % After reconstruction, the time axis is again of the form i \Delta t, with i>=1. Note that the first time bin is *before* measurement 1.
legend('True spikes', 'inferred spikes','binarized spikes');
set(gca,'FontSize',16);
linkaxes([h1,h2], 'x');


%% Same experiment, but without knowing the generative model parameters (blind).
Oalg.adaptive = 1;

tic;
[Ninf2,Cinf2,Palg2,Pphys2,Oalg2]=BSD( downF , Oalg);
toc;
binNinf2 = Ninf2 > Pphys2.threshold;



display(Pphys2)

figure(1); 
h1 = subplot(2,1,1);
plot(downT,downF,'Color',[0 0 0] + 0.5); hold on; 
plot(T,Cinf,'green'); plot(T,Cinf2,'blue'); legend('DF/F', 'Denoised DF/F','Denoised DF/F (blind)');
set(gca,'FontSize',16);
h2 = subplot(2,1,2);
plot(T,N,'black'); hold on;
plot(T,Ninf,'green'); 
plot(T,Ninf2,'blue');
legend('True spikes', 'inferred spikes','inferred spikes (blind)');
set(gca,'FontSize',16);
linkaxes([h1,h2], 'x');


%% Statistical Evaluation of the super-resolution: measure the response functions.

O.Time = 1e5;
T= [1:O.Time]/f;
[N,C,F] = BSD_generate_synthetic_signal(P,O);
SR = 5;
downF = F(SR:SR:end); % Down-sample the signal. 
downT = T(SR:SR:end); 

Oalg = struct; % Struct of experimental conditions & decoding options.
Oalg.Time = length(downF); % Number of time frames.
Oalg.dt = O.dt*SR; % interval duration. (s)
Oalg.nNeurons = O.nNeurons; % Number of neurons.
Oalg.adaptive = 0; % Not adaptive. Will use provided values for parameters, and estimate the unknown ones.
Oalg.superResolution = SR; % Attempt to reconstruct at the original signal frequency.

Palg=P; % Use known generative parameters

tic;
Ninf_SR =BSD( downF , Oalg , Palg); % Infer a spike signal with SR reconstructed at 100 Hz.
toc;

Oalg_no_SR = Oalg;
Oalg_no_SR.superResolution =1;

tic;
tmp = BSD( downF , Oalg_no_SR , Palg); % Infer a spike signal without SR at 20 Hz.
toc;

Ninf_no_SR = zeros(O.Time,1); % Transform to a 100 Hz signal by turning it into a piece-wise constant signal.

for i=1:O.Time;
    Ninf_no_SR(i) = tmp( ceil(i/SR) )/SR; % Such that tmp(i) = \sum_{j= (i-1) * SR + 1}^{i * SR} Ninf_no_SR(j)
end;


maxlag = 100;
PSF = estimate_PSF(N,Ninf_SR,maxlag); % infer point spread functions.
PSF_no_SR = estimate_PSF(N,Ninf_no_SR,maxlag); % infer point spread functions.
figure;
plot([-maxlag:maxlag] * O.dt, PSF,'LineWidth',2); hold on;
plot([-maxlag:maxlag] * O.dt, PSF_no_SR,'LineWidth',2);
xlim([-0.1,0.1]);
xlabel('Offset to spike (s)');
ylabel('Response function');
legend('Super-Resolution','Regular');
title(sprintf('Response function: f = %.f Hz, \\sigma = %.2f',f/SR,P.sigma));
set(gca,'FontSize',16);


%% Same, blind.

O.Time = 1e5;
T= [1:O.Time]/f;
[N,C,F] = BSD_generate_synthetic_signal(P,O);
SR = 5;
downF = F(SR:SR:end); % Down-sample the signal. 
downT = T(SR:SR:end); 

Oalg = struct; % Struct of experimental conditions & decoding options.
Oalg.Time = length(downF); % Number of time frames.
Oalg.dt = O.dt*SR; % interval duration. (s)
Oalg.nNeurons = O.nNeurons; % Number of neurons.
Oalg.adaptive = 1; % Not adaptive. Will use provided values for parameters, and estimate the unknown ones.
Oalg.superResolution = SR; % Attempt to reconstruct at the original signal frequency.

tic;
[Ninf_SR,~,~,Pphys] =BSD( downF , Oalg); % Infer a spike signal with SR reconstructed at 100 Hz.
toc;
display(Pphys);
Oalg_no_SR = Oalg;
Oalg_no_SR.superResolution =1;

tic;
tmp = BSD( downF , Oalg_no_SR); % Infer a spike signal without SR at 20 Hz.
toc;

Ninf_no_SR = zeros(O.Time,1); % Transform to a 100 Hz signal by turning it into a piece-wise constant signal.

for i=1:O.Time;
    Ninf_no_SR(i) = tmp( ceil(i/SR) )/SR; % Such that tmp(i) = \sum_{j= (i-1) * SR + 1}^{i * SR} Ninf_no_SR(j)
end;


maxlag = 100;
RF_blind = estimate_PSF(N,Ninf_SR,maxlag); % infer point spread functions.
RF_no_SR_blind = estimate_PSF(N,Ninf_no_SR,maxlag); % infer point spread functions.
figure;
plot([-maxlag:maxlag] * O.dt, PSF,'LineWidth',2); hold on;
plot([-maxlag:maxlag] * O.dt, PSF_no_SR,'LineWidth',2);
plot([-maxlag:maxlag] * O.dt, RF_blind,'LineWidth',2); hold on;
plot([-maxlag:maxlag] * O.dt, RF_no_SR_blind,'LineWidth',2);
xlim([-0.1,0.1]);
xlabel('Offset to spike (s)');
ylabel('Response function');
legend('Super-resolution','Regular','Super-Resolution (blind)','Regular (blind)');
title(sprintf('Response function (blind): f = %.f Hz, \\sigma = %.2f',f/SR,P.sigma));
set(gca,'FontSize',16);




