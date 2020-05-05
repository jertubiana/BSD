%% Figure 4b. 
% Comparison of the point-spread function width of inferred spike trains with and without super-resolution.

addpath(genpath('../'))
warning('off','MATLAB:nearlySingularMatrix')
P = struct; 
P.a = 1; P.b = 0; P.nu = 0.5; P.tauRaise = 0.1; P.tauDecay = 0.5;
O = struct;
fsignal = 5e2;
O.dt = 1/fsignal;
O.nNeurons= 1;
O.Time = 1e5; % 5 minutes recording.
f_acs = [1:fsignal];
f_acs = unique(fsignal./round(fsignal./f_acs));
sigmas = [0.01, 0.1, 0.2];
n_freq = length(f_acs);
n_sigmas = length(sigmas);
n_trials = 200;
maxlag = 1000;

PSF = zeros(n_freq,n_sigmas,n_trials,2*maxlag+1);
PSF_sr = zeros(n_freq,n_sigmas,n_trials,2*maxlag+1);
PSF_srb = zeros(n_freq,n_sigmas,n_trials,2*maxlag+1);
PSF_srbt = zeros(n_freq,n_sigmas,n_trials,2*maxlag+1);
tauRaises = zeros(n_freq,n_sigmas,n_trials);
tauDecays = zeros(n_freq,n_sigmas,n_trials);
tauRaisest = zeros(n_freq,n_sigmas,n_trials);
tauDecayst= zeros(n_freq,n_sigmas,n_trials);


for l = 1:n_trials
    for i = 1: n_sigmas
        P.sigma = sigmas(i);    
        [N,C,F] = BSD_generate_synthetic_signal(P,O);
        parfor j = 1: n_freq
            warning('off','MATLAB:nearlySingularMatrix')
            warning('off','MATLAB:SingularMatrix');
            f_ac = f_acs(j);    
            superResolution = round(fsignal/f_ac);
            downF = F(superResolution:superResolution:end);
            Oalg = O; 
            Oalg.slowOptimizer = 1;
            Oalg.dt = O.dt * superResolution;
            Oalg.Time = length(downF);
            Oalg.superResolution = superResolution;
            Oalg.adaptive = 0;
            fprintf('~~~~ Known params  %d %d %d~~~~\n',l,i,j);
            Ninf = BSD(downF,Oalg,P);
            fprintf('~~~~ No SR %d %d %d~~~~\n',l,i,j);
            Oalg.superResolution =1;
            tmp = BSD(downF,Oalg,P);
            
            Ninf2 = zeros(size(Ninf)); % Transform to a 100 Hz signal by turning it into a piece-wise constant signal.

            for k=1:O.Time;
                Ninf2(k) = tmp( ceil(k/superResolution) )/superResolution; % Such that tmp(i) = \sum_{j= (i-1) * SR + 1}^{i * SR} Ninf_no_SR(j)
            end;

            if length(Ninf)>length(N);
                Ninf = Ninf(1:length(N));
                Ninf2 = Ninf2(1:length(N));
            end;        
        
            PSF_sr(j,i,l,:) = estimate_PSF(N,Ninf,maxlag);
            PSF(j,i,l,:) = estimate_PSF(N,Ninf2,maxlag);        
        end;
    end;
end;

PSF_sr = squeeze(mean(PSF_sr,3));
PSF = squeeze(mean(PSF,3));


bell = @(params,X) params(1)*exp(-(X-params(3)).^2/params(2).^2);
options = optimoptions(@lsqcurvefit,'TolFun',1e-20,'Display','off','algorithm','levenberg-marquardt');


width_PSFsr = zeros(n_freq,n_sigmas);
widthPSF = zeros(n_freq,n_sigmas);

offset_PSFsr = zeros(n_freq,n_sigmas);
offset_PSF = zeros(n_freq,n_sigmas);

lag= [-maxlag:maxlag]'/fsignal;
for i = 1:n_sigmas
    parfor j = 1: n_freq
        x = lsqcurvefit(bell,[max(squeeze(PSF_sr(j,i,:))),O.dt,0],lag,squeeze(PSF_sr(j,i,:)),[],[],options);
        width_PSFsr(j,i) = abs(x(2));
        offset_PSFsr(j,i) = x(3);     

        x = lsqcurvefit(bell,[max(squeeze(PSF(j,i,:))),O.dt,0],lag,squeeze(PSF(j,i,:)),[],[],options);        
        widthPSF(j,i) = abs(x(2));
        offset_PSF(j,i) = x(3);        
    end;
end;
   
save('data_figure4b.mat');



%%

load('data_figure4b.mat');
Legend = cell(n_sigmas,1);
for i = 1:n_sigmas
    Legend{i} = sprintf('\\sigma = %.2f',sigmas(i));
end;

Legend{n_sigmas+1} = 'Regular';
Legend{n_sigmas+2}  = 'Super-Resolution';


figure('Position',[300,300,800,800]); plot(f_acs,1./width_PSF,'LineWidth',3); hold on;  
h1 = plot([0],[0],'Color','Black');
h2 = plot([0],[0],'--','Color','Black');


legend(Legend,'location','southeast')
set(gca,'ColorOrderIndex',1)
plot(f_acs,1./width_PSFsr,'--','LineWidth',3);
xlabel('f (Hz)');
ylabel('\delta t^{-1} (Hz)');
title('Response Function Inverse Width');
set(gca,'FontSize',28);
set(gcf, 'Color', 'w');
export_fig('figure4b.png','-nocrop');
export_fig('figure4b.eps','-nocrop');

