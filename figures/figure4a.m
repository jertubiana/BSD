clear all;
O = struct; % Struct of experimental conditions & decoding options.
O.Time = 1e2; % Number of time frames.
O.nNeurons =1; % Number of neurons.
f = 20; % Frequency of acquisition. (Hz)
O.dt = 1/f; % interval duration. (s)

P = struct; % Struct of generative model properties.
P.tauRise = 0.1; % Fluorescence raise time (s)
P.tauDecay = 0.5; % Fluorescence decay time (s)
P.nu = 0.5; % Neuron firing rate (Hz).
P.a = 1; % Spike amplitude
P.b = 0; % Baseline position
P.sigma = 1e-1; % Noise standard deviation.

[N,C,F] = BSD_generate_synthetic_signal(P,O);

superResolution = 4;
downF = F(superResolution:superResolution:end);

Oalg = struct;
Oalg.Time = length(downF);
Oalg.nNeurons = 1;
Oalg.dt = O.dt * superResolution;
Oalg.adaptive = 0;
Oalg.superResolution = superResolution;

Palg = struct;
Palg.tauRise = P.tauRise;
Palg.tauDecay = P.tauDecay;
Palg.sigma = P.sigma;
Palg.a = P.a;
Palg.b = P.b;
%%
tic;
[Ninf,Cinf,pAlg,pPhys,Oalg]=BSD( downF , Oalg , Palg); % Refer to BSD main code for detail of the output.
toc;

save('data_figure4a.mat');
%%
load('data_figure4a.mat');
T = [0:O.Time-1] * O.dt;

figure('Position',[300,300,800,800]); plot(T,F,'Color',[0 0 0] + 0.5,'LineWidth',1.5);  hold on;
plot(([superResolution:superResolution:Oalg.sTime]-1) * O.dt,downF,'LineWidth',2,'Color','blue');
plot(T,Cinf,'LineWidth',1.5,'Color','green');
plot(T(find(N)),-0.5*ones(size(find(N))),'square','Color','Black','MarkerSize',8,'MarkerFaceColor',[0 0 0]); 
plot(T,Ninf-1,'LineWidth',1.5,'Color','red')
xlabel('Time (s)');
title('Super-Resolution Deconvolution')
legend('Fluorescence','Sampled F','Reconstruction','spikes','inferred spikes');
set(gca,'FontSize',28);
set(gca,'Ytick',[]);
savefig('figure4a.fig');
set(gcf,'Color','White');
export_fig('figure4a.png','-nocrop');
export_fig('figure4a.eps','-nocrop');
