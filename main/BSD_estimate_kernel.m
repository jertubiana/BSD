function [tauRise,tauDecay] = BSD_estimate_kernel(F,sigma, O,varargin)


nsteps = O.nStepsKernelInference;
dt = O.dt;
if isfield(O,'tauRiseMin');
    tauRiseMin = 1.001 * O.tauRiseMin;
else
    tauRiseMin = 0;
end;

if isfield(O,'tauDecayMin');
    tauDecayMin = 1.001 * O.tauDecayMin;
else
    tauDecayMin = 0;
end;

if isfield(O,'tauRiseMax');
    tauRiseMax = 0.999 * O.tauRiseMax;
else
    tauRiseMax = inf;
end;

if isfield(O,'tauDecayMax');
    tauDecayMax= 0.999 * O.tauDecayMax;
else
    tauDecayMax = inf;
end;


autocov = xcov(F,F,nsteps,'biased');
sym_autocov = autocov(nsteps+1:end);
sigma2 = min(sigma^2, min(eig(toeplitz(sym_autocov)))); % make sure that sym_autocov is still a covariance matrix, ie positive definite.
sym_autocov(1) = sym_autocov(1) - sigma2;
sym_autocov = sym_autocov/sym_autocov(1);
sym_autocov(sym_autocov<0) = 0; % negative correlations come from either noise or regular spike trains with refractory periods; they are not accounted in the model.

options = optimoptions(@fmincon, 'display','off');
%options =  optimoptions(@fminunc,'display','off','Algorithm','quasi-newton');
if (O.est_tauRise == 1) && (O.est_tauDecay == 1);
    err = @(params) error_autocov(sym_autocov, params(1),params(2),dt,nsteps,O.superResolution);

    X0 = [0.1, 1];    
    A = zeros(1,2);
    A(1,1) = 1;
    A(1,2) = -1;
    B = -dt/10;
    [X,val,~] = fmincon(err,X0,A,B,[],[],[tauRiseMin,tauDecayMin],[tauRiseMax,tauDecayMax],[],options);
    %[X,val,~] = fminunc(err,X0,options);
    tauRise = min(X);
    tauDecay = max(X);
elseif (O.est_tauRise == 0) && (O.est_tauDecay == 1);
    tauRise = varargin{1};
    err = @(params) error_autocov(sym_autocov, tauRise,params(2),dt,nsteps,O.superResolution);
    X0 = 1;
    A = -1;
    B = -tauRise;
    [X,val,~] = fmincon(err,X0,A,B,[],[],[tauDecayMin],[tauDecayMax],[],options);
%     [X,val,~] = fminunc(err,X0,options);
    tauRise = min(tauRise,X);
    tauDecay = max(tauRise,X);
elseif (O.est_tauRise == 1) && (O.est_tauDecay == 0)
    tauDecay = varargin{1};
    err = @(params) error_autocov(sym_autocov, params(1),tauDecay,dt,nsteps,O.superResolution);
    X0 = min(0.1, tauDecay/2);    
    A = 1;
    B = tauDecay;
    [X,val,~] = fmincon(err,X0,A,B,[],[],[tauRiseMin],[tauRiseMax],[],options);
%     [X,val,~] = fminunc(err,X0,options);
    tauRise = min(X,tauDecay);
    tauDecay = max(X,tauDecay);
end
    
    


function error = error_autocov(sym_autocov, tauRise,tauDecay,dt,nsteps,superResolution)

lambdad = exp(-dt/tauDecay);
lambdar = exp(-dt/tauRise);
d = [0:nsteps]';
theoretical_autocov = lambdad.^d * psi(lambdad^2,superResolution)/(1-lambdad^2) - (lambdad.^d + lambdar.^d) * psi(lambdad * lambdar,superResolution)/(1-lambdad*lambdar) + lambdar.^d * psi(lambdar^2,superResolution)/(1-lambdar^2);
theoretical_autocov = theoretical_autocov/theoretical_autocov(1);
error = sum((sym_autocov - theoretical_autocov).^2);
end


function out = psi(lambda,superResolution)
    out = 1/superResolution * (lambda- 1)/(1-1/lambda^(1/superResolution));
end

end