%% Script for figure 1 of main article (without constrained-oopsi)



addpath(genpath('../'));

BSD_functions

O = struct; % Struct of experimental conditions & decoding options.
O.Time = 500; % Number of time frames.
O.nNeurons =1; % Number of neurons.
f = 10; % Frequency of acquisition. (Hz)
O.dt = 1/f; % interval duration. (s)

P = struct; % Struct of generative model properties.
P.tauRise = 0.1; % Fluorescence raise time (s)
P.tauDecay = 0.5; % Fluorescence decay time (s)
P.nu = 0.1; % Neuron firing rate (Hz).
P.a = 1; % Spike amplitude
P.b = 0; % Baseline position
P.sigma = 0.4; % Noise standard deviation.

tic;
[N,C,F] = BSD_generate_synthetic_signal(P,O);
toc;

Oalg = O;
Oalg.adaptive = 0;
Palg = struct; 
Palg.b = P.b;
Palg.tauRise = P.tauRise;
Palg.tauDecay = P.tauDecay;
tic;
[Ninf_BSD,Cinf_BSD,~,~] = BSD(F,Oalg,Palg);
toc;

Palg.lambda = P.sigma^2/(P.a * P.nu * O.dt);
tic;
[Ninf_oopsi,Cinf_oopsi,~,~] = BSD(F,Oalg,Palg);
toc;

P.gamma = gammaCoeff(P.tauRise,P.tauDecay,O.dt);
P.delta = deltaCoeff(P.tauRise,P.tauDecay,O.dt);
P.eta= etaCoeff(P.tauRise,P.tauDecay,O.dt);

% g = [P.gamma + P.delta, -P.delta];
% b = 0;
% c1 = [];
% tic;
% [~,~,~,~,~,Ninf_coopsi] = constrained_foopsi(F,b,c1,g);
% toc;

tic;
Ninf_naive = P.eta/P.a * filter([1 -P.gamma-P.delta,P.delta],1, F - P.b);
toc;


Oalg = O;
Oalg.adaptive = 0;
Palg = struct; 
Palg.b = P.b;
Palg.tauRise = P.tauRise;
Palg.tauDecay = P.tauDecay;
Palg.lambda = 0;
tic;
[Ninf_nothing,Cinf_nothing,~] = BSD(F,Oalg,Palg);
toc;

save('data_figure1.mat');




%% Load and plot

load('data_figure1.mat');




% legend('Fluorescence','non-negative','oopsi','con-oopsi','BSD')
legend('Fluorescence','non-negative','oopsi','BSD')
legend boxoff

Spikes = N;
Fluorescence = F;
% Calcium = {Cinf_nothing,Cinf_oopsi,Cinf_coopsi,Cinf_BSD};
% SpikesInferred = {Ninf_naive,Ninf_nothing,Ninf_oopsi,Ninf_coopsi,Ninf_BSD};
% Legends = {'naive','non-negative','oopsi','con-oopsi','BSD'};

Calcium = {Cinf_nothing,Cinf_oopsi,Cinf_BSD};
SpikesInferred = {Ninf_naive,Ninf_nothing,Ninf_oopsi,Ninf_BSD};
Legends = {'naive','non-negative','oopsi','BSD'};



T = [0:length(Fluorescence)-1] * O.dt;
lim = [T(1),T(end)];
T2 = T - lim(1);
lim2 = lim - lim(1);
figure('Position',[300,300,700,900]);
ax1=subplot(2,1,1);
plot(T2,Fluorescence,'Color',[0,0,0]+0.75,'LineWidth',2.0); hold on;
set(gca,'ColorOrderIndex',1)
for u = 1: length(Calcium);
    plot(T2,Calcium{u},'LineWidth',2.0);%,colors{u});
end;
legend('Fluorescence');
legend boxoff

plot(T2(find(Spikes)),0*ones(size(find(Spikes))),'square','Color','Black','MarkerSize',4,'MarkerFaceColor',[0 0 0]); hold on; 
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xlim(lim2);
ylabel('Original Signal');
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'FontSize',20)
box off
ax2=subplot(2,1,2);
hold on;
plot(T2(find(Spikes)),0*ones(size(find(Spikes))),'square','Color','Black','MarkerSize',4,'MarkerFaceColor',[0 0 0]); hold on; 
set(gca,'ColorOrderIndex',1)

nTrials = length(SpikesInferred);
for u = 1: length(SpikesInferred);
    if u ==1;
        plot(T2,SpikesInferred{u}+1.0+6*(nTrials-u)/(nTrials-1),'LineWidth',2.0,'Color',[0 0 0] + 0.75);
        set(gca,'ColorOrderIndex',1)
    else
        plot(T2,SpikesInferred{u}+6*(nTrials-u)/(nTrials-1),'LineWidth',2.0);
    end;
end;
legend('Spikes',Legends{:},'Location' ,[0.7,0.45,0.2,0.1]);
legend boxoff
xlim(lim2);
xlabel('Time (s)')
ylabel('Deconvolved Signal');
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'FontSize',20)

set(gcf, 'Color', 'w');
linkaxes([ax1,ax2],'x');
