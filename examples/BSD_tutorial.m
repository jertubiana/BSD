%% Basic test of the inference algorithm
% An introductory script to test BSD on synthetic fluorescence data.
% Execute each cell sequentially.

%% Generate a fluorescence signal

O = struct; % Struct of experimental conditions & decoding options.
O.Time = 1e3; % Number of time frames.
O.nNeurons =1; % Number of neurons.
f = 20; % Frequency of acquisition. (Hz)
O.dt = 1/f; % interval duration. (s)

P = struct; % Struct of generative model properties.
P.tauRise = 0.1; % Fluorescence raise time (s)
P.tauDecay = 0.5; % Fluorescence decay time (s)
P.nu = 0.5; % Neuron firing rate (Hz).
P.a = 1; % Spike amplitude
P.b = 0; % Baseline position
P.sigma = 0.2; % Noise standard deviation.


[N,C,F] = BSD_generate_synthetic_signal(P,O); % Generate synthetic signal


figure; plot(F,'Color',[0 0 0] + 0.5); hold on; plot(C); plot(N); legend('DF/F','Convolved Spikes','Spikes');
xlabel('Frames number');
title('Example of Synthetic signal');
set(gca,'FontSize',16)

%% Perform inference with known generative parameters.
Oalg = struct; % Struct of experimental conditions & decoding options.
Oalg.Time = O.Time; % Number of time frames.
Oalg.dt = O.dt; % interval duration. (s)
Oalg.nNeurons = O.nNeurons; % Number of neurons.
Oalg.adaptive = 0; % Not adaptive. Will use provided values for parameters, and estimate the unknown ones.

Palg=P; % Use known generative parameters

tic;
[Ninf,Cinf,pAlg,Pphys,Oalg]=BSD( F , Oalg , Palg); % Performs deconvolution.
toc;

binNinf = Ninf > Pphys.threshold;

figure(1); 
h1 = subplot(2,1,1);
plot(F,'Color',[0 0 0] + 0.5); hold on; plot(Cinf,'green'); legend('DF/F', 'Denoised DF/F');
set(gca,'FontSize',16);
h2 = subplot(2,1,2);
plot(N,'black'); hold on; plot(Ninf,'green'); plot(binNinf*0.5);
legend('True spikes', 'inferred spikes','binarized spikes');
set(gca,'FontSize',16);
linkaxes([h1,h2], 'x');

%% Now, perform inference without knowing the generative model parameters.

Oalg = struct; % Struct of experimental conditions & decoding options.
Oalg.Time = O.Time; % Number of time frames.
Oalg.dt = O.dt; % interval duration. (s)
Oalg.nNeurons = O.nNeurons; % Number of neurons.
Oalg.adaptive = 0; % Not adaptive. Estimate the generative parameters once, then deconvolve.

tic; 
[Ninf,Cinf,Palg,Pphys,Oalg]=BSD( F , Oalg);
toc;

display(Pphys)

binNinf = Ninf > Pphys.threshold;

figure(2); 
h1 = subplot(2,1,1);
plot(F,'Color',[0 0 0] + 0.5); hold on; plot(Cinf,'green'); legend('DF/F', 'Denoised DF/F');
set(gca,'FontSize',16);
h2 = subplot(2,1,2);
plot(N,'black'); hold on; plot(Ninf,'green'); plot(binNinf*0.5);
legend('True spikes', 'inferred spikes','binarized spikes');
set(gca,'FontSize',16);
linkaxes([h1,h2], 'x');


%% Some times 1 iteration is not enough. Refine iteratively the parameters.

Oalg = struct; % Struct of experimental conditions & decoding options.
Oalg.Time = O.Time; % Number of time frames.
Oalg.dt = O.dt; % interval duration. (s)
Oalg.nNeurons = O.nNeurons; % Number of neurons.
Oalg.adaptive = 1; % Adaptive. Will use the provided values in P as initializer, estimate the unprovided ones, then iteratively refine the values.
Oalg.iterations = 10; % Maximal number of iterations. Default: 5.

tic; 
[Ninf2,Cinf2,Palg,Pphys,Oalg]=BSD( F , Oalg);
toc;

display(Pphys)
close all;

figure(3); 
h1 = subplot(2,1,1);
plot(F,'Color',[0 0 0] + 0.5); hold on; plot(Cinf,'red'); plot(Cinf2,'green'); legend('DF/F', 'Denoised DF/F 1 iteration',sprintf('Denoised DF/F %d iterations',Oalg.iterations));
set(gca,'FontSize',16);
h2 = subplot(2,1,2);
plot(N,'black'); hold on; plot(Ninf,'red'); plot(Ninf2,'green'); legend('Spikes', 'Inferred Spikes 1 iteration',sprintf('Inferred Spikes %d iterations',Oalg.iterations));
set(gca,'FontSize',16);
linkaxes([h1,h2], 'x');



%% ~~~~~ Typical use case: Suppose the raise and decay time scales are known.

O = struct; % Struct of experimental conditions & decoding options.
O.Time = 1e3; % Number of time frames.
O.nNeurons =1; % Number of neurons.
f = 20; % Frequency of acquisition. (Hz)
O.dt = 1/f; % interval duration. (s)

P = struct; % Struct of generative model properties.
P.tauRise = 0.1; % Fluorescence raise time (s)
P.tauDecay = 0.5; % Fluorescence decay time (s)
P.nu = 0.5; % Neuron firing rate (Hz).
P.a = 1; % Spike amplitude
P.b = 0; % Baseline position
P.sigma = 0.2; % Noise standard deviation.
[N,C,F] = BSD_generate_synthetic_signal(P,O);

Oalg = O;
Oalg.adaptive = 1; % Adaptive; => iterate over b
Oalg.est_tauRise = 0; % Stick with the default tauRise
Oalg.est_tauDecay = 0; % Stick with the default tauDecay

% ~~ Other example with other variables ~~ 
% Palg.lambda = my_lambda % Use your own sparsity prior
% Oalg.est_lambda = 0;

% Palg.b = my_baseline; % Use your own baseline.
% Oalg.est_b = 0;

% Palg.sigma = my_sigma; % Use your own noise level.
% Oalg.est_sigma = 0; 

% ~~~~~~~~~~

Palg=struct;
Palg.tauRise = 0.1;
Palg.tauDecay = 0.5;
[Ninf,Cinf,Palg,Pphys,Oalg]=BSD( F , Oalg , Palg);

binNinf = Ninf > Pphys.threshold;

figure(2); 
h1 = subplot(2,1,1);
plot(F,'Color',[0 0 0] + 0.5); hold on; plot(Cinf,'green'); legend('DF/F', 'Denoised DF/F');
set(gca,'FontSize',16);
h2 = subplot(2,1,2);
plot(N,'black'); hold on; plot(Ninf,'green'); plot(binNinf*0.5);
legend('True spikes', 'inferred spikes','binarized spikes');
set(gca,'FontSize',16);
linkaxes([h1,h2], 'x');

%% ~~~~~ Typical use case: use parallel computing when multiple neurons. ~~~

O.nNeurons = 10;
[N,C,F] = BSD_generate_synthetic_signal(P,O);

Oalg = O;
[Ninf,Cinf,Palg,Pphys,Oalg]=pBSD( F , Oalg);

figure; plot(F); hold on; plot(C); plot(N); legend('Fluorescence','denoised Fluorescence','Spikes');

display(Pphys);