addpath(genpath('../'))
O = struct; % Struct of experimental conditions & decoding options.
O.Time = 5e3; % Number of time frames.
O.nNeurons =1; % Number of neurons.
f = 60; % Frequency of acquisition. (Hz)
O.dt = 1/f; % interval duration. (s)

P = struct; % Struct of generative model properties.
P.tauRise = 0.1; % Fluorescence raise time (s)
P.tauDecay = 0.5; % Fluorescence decay time (s)
P.nu = 0.5; % Neuron firing rate (Hz).
P.a = 1; % Spike amplitude
P.b = 0; % Baseline position
P.sigma = 0.5; % Noise standard deviation.


[N,C,F] = BSD_generate_synthetic_signal(P,O); % Generate synthetic signal


%% Perform inference with known generative parameters.
Oalg = struct; % Struct of experimental conditions & decoding options.
Oalg.Time = O.Time; % Number of time frames.
Oalg.dt = O.dt; % interval duration. (s)
Oalg.nNeurons = O.nNeurons; % Number of neurons.
Oalg.adaptive = 1; % Not adaptive. Will use provided values for parameters, and estimate the unknown ones.
Oalg.tauRiseMin = 0.05; % Inferring accuraate rise time can be challenging at high noise level. Provide a little help...

tic;
[Ninf,Cinf,pAlg,Pphys,Oalg]=BSD( F , Oalg); % Performs deconvolution.
toc;

[Ninfconvolved , PSF,DeltaT]= BSD_convolve_with_PSF(Ninf,Oalg,pAlg);


display('Zoom-in on each spike to see the effect of convolution!');
figure(1); 
h1 = subplot(2,1,1);
plot(F,'Color',[0 0 0] + 0.5); hold on; plot(Cinf,'green'); legend('DF/F', 'Denoised DF/F');
set(gca,'FontSize',16);
h2 = subplot(2,1,2);
plot(N,'black'); hold on; plot(Ninf,'green');
plot(Ninfconvolved,'red');
legend('True spikes', 'inferred spikes','inferred spikes (convolved)');
set(gca,'FontSize',16);
linkaxes([h1,h2], 'x');



figure(2);
hold on;
plot(DeltaT,PSF/max(PSF),'LineWidth',3);
for n = 0:5;
    plot([O.dt*n,O.dt*n],[0,1],'Color',[0.2,0.2,0.2]);
    plot([-O.dt*n,-O.dt*n],[0,1],'Color',[0.2,0.2,0.2]);
end;
xlabel('Time');
ylabel('PSF');
title('Point-spread function');
set(gca,'FontSize',16);