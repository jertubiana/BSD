
%% Script for generating figure 6 panels b,c,d.
addpath(genpath('../'))

q = 0.15;
tauBaselines = [60, 10, 60, 30, 40, 40, 10, 60, 30, 20]; % Values determined on train set.

index =1 ; % Panel choice.
if index == 1
    dataset = 1;
    neuron = 2;        
elseif index==2
    dataset = 1;
    neuron = 2;    
elseif index ==3
    dataset = 7;
    neuron = 1;
else
    dataset = 7;
    neuron = 1;    
end;



nAlgorithms = 4;
convolved = [false,false,false,true];

folders = {'../SpikeFinder/predictions_noregularization/',...
    '../SpikeFinder/predictions_notadaptive/',...
    '../SpikeFinder/predictions_groundtruth/',...
    '../SpikeFinder/predictions_groundtruth/'...    
    };

Legends__ = {'Non-negative','BSD','Adaptive BSD','Adaptive BSD * PSF'};

colors = {'red','green','blue','magenta'};

Calcium = cell(nAlgorithms,1);
SpikesInferred = cell(nAlgorithms,1);


%
load(sprintf('Datasets/%d.train.mat',dataset) );
Fluorescence = normalize_remove_baseline(F{neuron},q,tauBaselines(dataset),O.dt);
Spikes =  N{neuron};


for u = 1:nAlgorithms;
    folder = folders(u);
    is_convolved = convolved(u);
    load(strcat(folder,sprintf('%d.train.mat',dataset) ) )
    if is_convolved
        SpikesInferred{u} = Ninf_convolved{neuron};
    else
        SpikesInferred{u} = Ninf{neuron};
    end;
    Calcium{u} = Cinf{neuron};
end;


%% Display results.

T = [0:length(Fluorescence)-1] * dt(neuron);

lim = [T(1),T(end)];
if index == 1
    lim = [321,323]; % Dataset 1, neuron 2.
    lim_zoom =  [321.6,321.7];
    position = [.285 .42 .2 .15]; % num 1. 
    hbracket = 8.5;
elseif index == 2
    lim = [149,153]; % Dataset 1, neuron 2.
    lim_zoom = [150.75,150.88]; % Dataset 1, neuron 2.
    position = [.375 .41 .2 .17]; % num 3.
    hbracket = 10.25;
    
elseif index ==3
    lim = [104,107]; % Dataset 7, neuron 1.
    lim_zoom = [105.02,105.1];
    position = [.30 .46 .2 .14]; % num 3.
    hbracket = 10;
end;

T2 = T - lim(1);
lim2 = lim - lim(1);
lim_zoom2 = lim_zoom - lim(1);
figure('Position',[300,300,700,900]);
ax1=subplot(2,1,1);
plot(T2,Fluorescence,'Color',[0,0,0]+0.25,'LineWidth',2.0); hold on;
set(gca,'ColorOrderIndex',1)
for u = 1: nAlgorithms
    plot(T2,Calcium{u},'LineWidth',2.0);%,colors{u});
end;
plot(T2(find(Spikes)),0*ones(size(find(Spikes))),'square','Color','Black','MarkerSize',4,'MarkerFaceColor',[0 0 0]); hold on; 
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xlim(lim2);
ylabel('Original Signal');
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'FontSize',24)
box off
ax2=subplot(2,1,2);
hold on;
plot([-1,0],[0,0],'Color',[0,0,0]+0.25,'LineWidth',2.0);
plot(T2(find(Spikes)),0*ones(size(find(Spikes))),'square','Color','Black','MarkerSize',4,'MarkerFaceColor',[0 0 0]); hold on; 
set(gca,'ColorOrderIndex',1)
for u = 1: nAlgorithms;
    plot(T2,SpikesInferred{u}+6*(nAlgorithms-u)/(nAlgorithms-1),'LineWidth',2.0);
end;


plot([lim_zoom2(1),lim_zoom2(1)],[hbracket-0.5,hbracket ],'Color','Black','LineWidth',2.0);
plot([lim_zoom2(2),lim_zoom2(2)],[hbracket-0.5,hbracket ],'Color','Black','LineWidth',2.0);
plot([lim_zoom2(1),lim_zoom2(2)],[hbracket ,hbracket ],'Color','Black','LineWidth',2.0);
legend('Fluorescence','Spikes',Legends__{:},'Location' ,[0.7,0.4,0.2,0.1]);
legend boxoff
xlim(lim2);
xlabel('Time (s)')
ylabel('Deconvolved Signal');
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'FontSize',24)
set(gcf, 'Color', 'w');
linkaxes([ax1,ax2],'x');


axes('Position',position)
box on
hold on;
plot([-1,0],[0,0],'Color',[0,0,0]+0.25,'LineWidth',2.0);
all_locations = find(Spikes);
[~,I]=min(abs(T(all_locations) - mean(lim_zoom) ));
location = T2(all_locations(I));
plot([location, location],[0,5],'Color','black','LineWidth',2.0);
hold on;
plot(T2(find(Spikes)),0*ones(size(find(Spikes))),'square','Color','Black','MarkerSize',4,'MarkerFaceColor',[0 0 0]); hold on; 
set(gca,'ColorOrderIndex',1)
window_ = [all_locations(I)-10:all_locations(I)+10];
ticks = T2(window_);
tick_labels = cell(length(window_),1);
for u=1:length(window_);
    tick_labels{u} ='';
end;

for u = 1: nAlgorithms;
    maxi = max( SpikesInferred{u}(window_) );
    h=plot(T2,SpikesInferred{u}/maxi+4*(nAlgorithms-u)/(nAlgorithms-1),'LineWidth',2.0,'Marker','s');
    set(h, 'markerfacecolor', get(h, 'color'));
end;
xlim(lim_zoom2)
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'xtick',ticks)
set(gca,'xticklabel', tick_labels);
grid on


savefig(sprintf('figure6_%d.fig', index ) );
export_fig(sprintf('figure6_%d.png',index),'-nocrop');
