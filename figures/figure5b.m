%% Script for figure 5b: 
addpath(genpath('../'))
load('../SpikeFinder/ground_truth_time_constants.mat');



datasets = [1,6,7,8,9,10];

tauRises = tauRiseGroundTruth(datasets);
tauDecays = tauDecayGroundTruth(datasets);

names = {'OGB-1','GCaMP5k','GCaMP6f','GCaMP6s','jRCAMP1a','jRGECO1a'};

O = struct;
O.dt = 1/60;
O.z1 = 2.326;
O.z2 = 2.326;
O.delta_max = 20;
O.discretisation = 40;
O.Nsim = 1e5;

sigmas = 1e-4+ [0:0.025:1];


PSF = zeros(length(sigmas),length(tauRises),2*O.discretisation*O.delta_max+1);
width = zeros(length(sigmas),length(tauRises));
lag = zeros(length(sigmas),length(tauRises));
for k = 1:length(sigmas);
    display(k);
    for l = 1: length(tauRises);        
        P = struct;
        P.sigma = sigmas(k);
        P.tauRise = tauRises(l);
        P.tauDecay = tauDecays(l);
        P.a = 1;
        P.b = 0;
        [RF,DeltaT,lag(k,l),width(k,l),~,~] = BSD_theoretical_accuracy(P,O);
        PSF(k,l,:) = RF';
    end;
end;

save('data_figure5b.mat');
%%
load('data_figure5b.mat');


figure;
all_default_colors = get(gca,'colororder');
close;


colors = [4,2,1,3,6,7];


k0 = 10;
window = [801-100:801+100];

figure1=figure('Position',[300,300,800,800]);
hold on;
for l = 1:length(tauRises);
    plot(sigmas, width(:,l), 'LineWidth',3.,'Color',all_default_colors(colors(l),:));
end;
set(gca,'FontSize',18);
xlabel('$\frac{\sigma}{a}$','interpreter','latex');
ylabel('PSF width \delta t (s)');
xlim([0,1]);
title(sprintf('Width of point-spread function f = %.f Hz',1/O.dt))
legend(names{:},'Location','southeast');
legend boxoff
set(gca,'FontSize',26)
axes('Position',[.25 .6 .3 .3])
box on
hold on;
for l = 1:length(tauRises);
    plot(DeltaT(window),squeeze(PSF(k0,l,window)),'LineWidth',2,'Color',all_default_colors(colors(l),:));
end;
xlabel('\Delta t (s)');
ylabel('Lag probability');
set(gca,'ytick',[])
set(gca,'yticklabel',[])

% Create doublearrow
annotation(figure1,'doublearrow',[0.38 0.46],...
    [0.660375661375661 0.66005291005291]);

% Create textbox
annotation(figure1,'textbox',...
    [0.396 0.620370370370369 0.0565000000000001 0.0277777777777777],...
    'String','\delta t',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Helvetica',...
    'FitBoxToText','off');
set(gca,'FontSize',22);
set(gca,'YTick',[0,0.01,0.02,0.03,0.04])
set(gcf, 'Color', 'w');
savefig('figure_5b.fig');
export_fig('figure_5b.png','-nocrop');