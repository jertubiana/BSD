function [a,b,sigma, tauRise, tauDecay, gamma, delta, eta, lambda, mu, threshold]=BSD_initialization(Fluorescence,O,varargin)

% A function to get initial estimates of the spike trains and the
% parameters of the generative model of the problem. 
% Input: Fluorescence (Time X Neurons).
%        O (structure of algorithms options).
%       P (structure of algorithms parameters).
% Output: N Rough spike train estimate.
%         C: Corresponding convolved signal.
         % Pinit: Completed P.

         
%% Define some global functions.    
BSD_functions

%% Parse arguments
if length(varargin) == 0;
    P = struct;
elseif length(varargin) ==1;
    P = varargin{1};
    k = 1;
else
    P = varargin{1}; % structure of known parameters
    k = varargin{2}; % neuron index
end;


%% Initialization.

valid_index = ~isnan(Fluorescence);

if ~isfield(P,'b');
    nbin = max(10, round(length(Fluorescence(valid_index))/50));
    [alpha,beta]=hist(Fluorescence(valid_index),nbin); 
    b= beta(find(alpha==max(alpha),1,'first'));
else
    if length(P.b)>1
        b = P.b(k);
    else
        b = P.b;
    end;
end;

signal = Fluorescence - b;
signal(~valid_index) = 0;

if ~isfield(P,'sigma');
    tmp2 = [signal(signal<0); -signal(signal<0)]; % Fit the negative part of the signal as a gaussian distribution.
    sigma = std(tmp2);
    sigma = max(sigma,O.conditioner*100);
else 
    if length(P.b)>1
        sigma = P.sigma(k);
    else
        sigma = P.sigma;
    end;
end;
    

if (~isfield(P,'tauRise') ) || (~isfield(P,'tauDecay'));
    if isfield(P,'tauRise');
        if length(P.tauRise)>1
            tauRise = P.tauRise(k);
        else
            tauRise = P.tauRise;
        end;
        O.est_tauRise = 0; O.est_tauDecay =1;
        [tauRise, tauDecay] = BSD_estimate_kernel( signal,sigma, O, tauRise);
    elseif isfield(P,'Decay');
        if length(P.tauDecay)>1
            tauDecay= P.tauDecay(k);
        else
            tauDecay = P.tauDecay;
        end;
        O.est_tauDecay = 0; O.est_tauRise = 1;      
        [tauRise,tauDecay] = BSD_estimate_kernel( signal,sigma, O, tauDecay);
    else
        O.est_tauRise = 1; O.est_tauDecay = 1;
        [tauRise,tauDecay] = BSD_estimate_kernel( signal,sigma, O);
    end;
else
    if length(P.tauRise)>1
        tauRise = P.tauRise(k);
    else
        tauRise = P.tauRise;
    end;
    if length(P.tauDecay)>1
        tauDecay = P.tauDecay(k);
    else
        tauDecay = P.tauDecay;
    end;
end;
    
eta = etaCoeff(tauRise, tauDecay, O.dt/O.superResolution);
gamma = gammaCoeff(tauRise, tauDecay, O.dt/O.superResolution);
delta = deltaCoeff(tauRise, tauDecay, O.dt/O.superResolution);


if ~isfield(P,'a');
    tmp1 = (var(Fluorescence)-sigma^2); % ~ a^2 ||K||_2^2 nu dt
    tmp2 = (mean(Fluorescence-b)); % ~ a * nu * dt * ||K||_1
    if (tmp1>0) && (tmp2>0);
        if O.superResolution>1
            factor = mean(normKernelCoeffSR1(tauRise,tauDecay,O.dt,O.superResolution)./normKernelCoeffSR(tauRise,tauDecay,O.dt,O.superResolution).^2 ) ;
        else
            factor = normKernelCoeff1(tauRise,tauDecay,O.dt)/normKernelCoeff(tauRise,tauDecay,O.dt)^2;
        end;
        a = tmp1/tmp2 * factor;
        flag = 1;
    else
        a = 0;
        flag = 0;
    end;
    if (a< sigma/normKernelCoeff(tauRise,tauDecay,O.dt) ) || (a>5*max(Fluorescence-b)) || (flag ==0) % if SNR_eff < 1, theoretically impossible to find spikes. set SNR_eff to 1.
        a = sigma/max(1,normKernelCoeff(tauRise,tauDecay,O.dt) ); % a < sigma.
        flag = 0; % use the formula without a.
    end;
else
    flag = 1;
    if length(P.a)>1
        a = P.a(k);
    else
        a = P.a;
    end;
    
end;

    

if ~isfield(P,'lambda');
    if flag == 0; % estimate of a is unreliable; likely due to signal artifacts or very low SNR.
        if O.superResolution>1;
            lambda=sparsityCoeffSR(tauRise,tauDecay,O.dt,sigma, O.z1,O.superResolution);
        else
            lambda=sparsityCoeff(tauRise,tauDecay,O.dt,sigma, O.z1);
        end;
    else
        if O.superResolution>1;
            lambda=sparsityCoeffSR2(tauRise,tauDecay,O.dt,sigma, a,O.z1,O.z2,O.superResolution);
        else
            lambda=sparsityCoeff2(tauRise,tauDecay,O.dt,sigma,a, O.z1,O.z2);
        end;
    end;
else
    if length(P.lambda)>1
        lambda=P.lambda(k);
    else
        lambda = P.lambda;
    end;
end;

mu = sparsityCoeffSR(tauRise,tauDecay,O.dt,1,1,O.superResolution).^2; % The conditioner for the deconvolution. Must be proportional to |K_j| in order to have flat solution when tauR is small.


if ~isfield(P,'threshold'); 
    threshold=thresholdCoeff(tauRise, tauDecay,O.dt,sigma, O.z3);
else
    if length(P.threshold)>1
        threshold = P.threshold(k);
    else
        threshold = P.threshold;
    end;
end;

    
end