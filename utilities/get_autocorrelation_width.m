function width = get_autocorrelationwidth(Spike)

if size(Spike,1) > 1
    Spike = Spike';
end;

nmax = 200;

autocorrelation = xcorr(Spike,Spike,nmax,'unbiased') - mean(Spike)^2;
autocorrelation = autocorrelation/autocorrelation(nmax+1);


%% Solution 1.
% autocorrelation = 0.5* (autocorrelation(nmax+1:end) + autocorrelation(nmax+1:-1:1));
% width = sqrt( sum([0:nmax].^2 .* abs(autocorrelation) )/sum(abs(autocorrelation)) );
% figure; plot([0:nmax],autocorrelation); %hold on; plot([-nmax:nmax],bell(params_opt,[-nmax:nmax]));



%% Solution 2.
% options = optimoptions(@lsqcurvefit,'Display','off');
% bell = @(params,X) params(1)*exp(-X.^2/(2*params(3)^2)) + params(2);
% 
% maxi0 = 1;
% offset0 = 0.5 * ( mean(autocorrelation(1: nmax/10) ) + mean(autocorrelation(end-nmax/10+1: end) ) );
% % width0 = sqrt( sum([-nmax:nmax].^2 .* abs(autocorrelation) )/sum(abs(autocorrelation)) );
% width0 = 1;
% 
% params_opt = lsqcurvefit(bell,[maxi0,offset0,width0],[-nmax:nmax],autocorrelation,[],[],options);
% width = abs(params_opt(3));
% 
% figure; plot([-nmax:nmax],autocorrelation); hold on; plot([-nmax:nmax],bell(params_opt,[-nmax:nmax]));

%% Solution 3.
autocorrelation_mod = autocorrelation;
autocorrelation_mod(nmax+1) = 0.5 * (autocorrelation_mod(nmax) + autocorrelation_mod(nmax+2) );
options = optimoptions(@lsqcurvefit,'Display','off');
bell = @(params,X) params(1)*exp(-X.^2/(2*params(2)^2));

maxi0 = autocorrelation_mod(nmax+1);
width0 = 1;

params_opt = lsqcurvefit(bell,[maxi0,width0],[-nmax:nmax],autocorrelation_mod,[],[],options);
width = sqrt( params_opt(1)) * abs(params_opt(2));

% figure; plot([-nmax:nmax],autocorrelation); hold on; plot([-nmax:nmax],bell(params_opt,[-nmax:nmax]));


