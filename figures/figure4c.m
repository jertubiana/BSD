clear all
addpath(genpath('../'))

P = struct; 
P.a = 1; P.b = 0; P.nu = 1; P.tauRise = 0.1; P.tauDecay = 0.5;
O = struct;
fsignal = 5e2;
O.dt = 1/fsignal;
O.nNeurons= 1;
O.Time = 1e5;
Legend = cell(3,1);
sigmas = [0.01,0.05,0.2];
maxlag = 100;
ntrials = 50;
PSFs = zeros(ntrials,3,2*maxlag+1);
for i =1:3
    display(i)
    P.sigma = sigmas(i);
    parfor j=1:ntrials
        display(j)
        [N,C,F] = BSD_generate_synthetic_signal(P,O);      
        f_ac = 10;
        superResolution = fsignal/f_ac;        
        downF = F(superResolution:superResolution:end);          
        Oalg = O;
        Oalg.slowOptimizer = 1;
        Oalg.dt = O.dt * superResolution;
        Oalg.Time = length(downF);
        Oalg.superResolution = superResolution;
        Oalg.adaptive = 0;
        Ninf = BSD(downF,Oalg,P);
        PSFs(j,i,:)  = estimate_PSF(N,Ninf,maxlag);
    end;
end;
mPSF = squeeze(mean(PSFs,1));
save('data_figure4c.mat');
%%
load('data_figure4c.mat');
figure('Position',[300,300,800,800]);
hold on;plot([-maxlag:maxlag]/fsignal,mPSF,'LineWidth',3);
for i = 1:3
    Legend{i} = sprintf('\\sigma = %.2f',sigmas(i));
end;
legend(Legend{:},'FontSize',28);
xlabel('$\delta t$ (s)','interpreter','latex');
xlim([-0.1,0.1]);
ylim([-0.01,max(mPSF(:))])
title('Response function, f = 10 Hz');
set(gca,'FontSize',28);
savefig('figure4c.fig');
set(gcf, 'Color', 'w');
export_fig('figure4c.png','-nocrop');


