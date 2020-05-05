function [C,cost,bNew,aNew,sigmaNew,lambdaNew,muNew,thresholdNew,gammaNew,deltaNew,etaNew,tauRiseNew,tauDecayNew] ...
                 = BSD_parameter_estimation(N,Fluorescence,O,tauRise,tauDecay,b,a,sigma,threshold,lambda)    % update parameters based on previous iteration.     


% Function that refines the generative model parameters given the
% fluorescence signal, the spike train and the old parameters.
% Input: N: an individual Spike train.
%        Fluorescence: an individual fluorescence trace. 
%        O: the algorithm options.
%        tauRise,tauDecay,b,a,sigma,threshold, lambda: the old parameters.
% Output: 
%         C: the convolved signal corresponding to N C = K* N+b.
%         cost: the cost function value after parameter refinement.
%         bNew,
%         aNew,sigmaNew,lambdaNew,muNew,thresholdNew,gammaNew,deltaNew,etaNew,tauRiseNew,tauDecayNew:
%         the new generative model parameters.


%% Initialization.
BSD_functions

if isfield(O,'tauRiseMin');
    tauRiseMin = O.tauRiseMin;
else
    tauRiseMin = 0;
end;

if isfield(O,'tauDecayMin');
    tauDecayMin = O.tauDecayMin;
else
    tauDecayMin = 0;
end;

if isfield(O,'tauRiseMax');
    tauRiseMax = O.tauRiseMax;
else
    tauRiseMax = inf;
end;

if isfield(O,'tauDecayMax');
    tauDecayMax= O.tauDecayMax;
else
    tauDecayMax = inf;
end;

   


%% Refine the convolution kernel of the signal.
if O.est_tauDecay || O.est_tauRise
    nmax = max(min(round ( 10* (tauRise+tauDecay)/O.dt), floor(O.Time * 0.75) ),5);
    sF2 = sum(Fluorescence.^2);
    sF = sum(Fluorescence);
    if O.thresholdBeforeKernelInference ==1;
        N(N<=threshold) = 0;
    end;
        
    
    if O.superResolution >1
        sN = zeros(O.superResolution,1);
        for i =1:O.superResolution;
            sN(i) = sum(N(i:O.superResolution:end));
        end;
        
        FxN = zeros(nmax+1,O.superResolution);
        BxN = zeros(nmax+1,O.superResolution);
        for i =1:O.superResolution;
            tmp = xcorr(Fluorescence,N(i:O.superResolution:end),nmax);
            FxN(:,i) = tmp(nmax+1:end);
            tmp = xcorr(ones(size(Fluorescence)),N(i:O.superResolution:end),nmax);
            BxN(:,i) = tmp(nmax+1:end);
        end;        
        
        NxN = zeros(2*nmax+1,O.superResolution,O.superResolution);
        for i=1:O.superResolution;
            for j=1:O.superResolution;
                NxN(:,i,j) = xcorr(N(i:O.superResolution:end),N(j:O.superResolution:end),nmax);
            end;
        end;

        lam = lambda'./normKernelCoeffSR(tauRise,tauDecay,O.dt,O.superResolution);
        cost_function = @(params) BSD_fcostSR(params(1),params(2), params(3),sF, sF2,sN, BxN, FxN, NxN,N, lam, O.dt,O.Time,nmax,O.superResolution);
    else
        sN = sum(N);
        FxN = xcorr(Fluorescence,N, nmax);
        FxN = FxN(nmax+1:end);
        BxN = xcorr(ones(size(N)),N, nmax);
        BxN = BxN(nmax+1:end);        
        NxN = xcorr(N,N,nmax);
        lam = lambda/normKernelCoeff(tauRise,tauDecay,O.dt);
        cost_function = @(params) BSD_fcost(params(1),params(2), params(3), sF, sF2,sN, BxN,FxN, NxN,N,lam,O.dt,O.Time,nmax);
    end;
    
%     options =  optimoptions(@fminunc,'display','off','Algorithm','quasi-newton');    
    options = optimoptions(@fmincon, 'display','off');
    if (O.est_tauDecay) && (O.est_tauRise) && (O.est_b);
        X0 = [tauRise,tauDecay,b];
        A = zeros(1,3);
        A(1,1) = 1;
        A(1,2) = -1;
        B = -O.dt/10;         
        [X,cost,~] = fmincon(cost_function,X0,A,B,[],[],[tauRiseMin,tauDecayMin,-inf],[tauRiseMax,tauDecayMax,+inf],[],options);                 
%         [X,cost] = fminunc(cost_function,X0,options);
        tauRiseNew = min(X(1:2));
        tauDecayNew = max(X(1:2));
        bNew = X(3);
    elseif (O.est_tauDecay) && (O.est_tauRise) && (~O.est_b);
        cost_function = @(params) cost_function([params(1),params(2),b]);
        X0 = [tauRise,tauDecay];
        [X,cost,~] = fminunc(cost_function,X0,options);
        tauRiseNew = min(X);
        tauDecayNew = max(X);
        bNew = b;
    elseif (O.est_tauDecay) && (~O.est_tauRise) && (O.est_b);
        cost_function = @(params) cost_function([tauRise,params(1),params(2)]);
        X0 = [tauDecay,b];
        [X,cost,~] = fminunc(cost_function,X0,options);
        tauRiseNew = tauRise;
        tauDecayNew = X(1);
        bNew = X(2);
    elseif (O.est_tauDecay) && (~O.est_tauRise) && (~O.est_b);
        cost_function = @(params) cost_function([tauRise,params,b]);
        X0 = tauDecay;
        [X,cost,~] = fminunc(cost_function,X0,options);
        tauRiseNew = tauRise;
        tauDecayNew = X;
        bNew = b;
    elseif (~O.est_tauDecay) && (O.est_tauRise) && (O.est_b);
        cost_function = @(params) cost_function([params(1),tauDecay,params(2)]);
        X0 = [tauRise,b];
        [X,cost,~] = fminunc(cost_function,X0,options);
        tauRiseNew = X(1);
        tauDecayNew = tauDecay;
        bNew = X(2);
    elseif (~O.est_tauDecay) && (O.est_tauRise) && (~O.est_b);
        cost_function = @(params) cost_function([params,tauDecay,b]);
        X0 = tauRise;
        [X,cost,~] = fminunc(cost_function,X0,options);
        tauRiseNew = X;
        tauDecayNew = tauDecay;
        bNew = b;        
    end;
    tauRiseNew = real(tauRiseNew); % just in case.
    tauDecayNew = real(tauDecayNew);
%     tauRiseNew = max(min(tauRiseNew,O.tauRiseMax),O.tauRiseMin);
%     tauDecayNew = max(min(tauDecayNew,O.tauDecayMax),O.tauDecayMin);

    % Update the discrete dynamics parameters
    deltaNew = deltaCoeff(tauRiseNew,tauDecayNew,O.dt/O.superResolution);
    gammaNew = gammaCoeff(tauRiseNew,tauDecayNew,O.dt/O.superResolution);
    etaNew = etaCoeff(tauRiseNew,tauDecayNew,O.dt/O.superResolution);
    C=filter(1,[etaNew -etaNew*(gammaNew+deltaNew),etaNew*deltaNew],N);  % Update C.
    C = C + bNew;
else
    tauRiseNew = tauRise;
    tauDecayNew = tauDecay;
    deltaNew = deltaCoeff(tauRiseNew,tauDecayNew,O.dt/O.superResolution);
    gammaNew = gammaCoeff(tauRiseNew,tauDecayNew,O.dt/O.superResolution);
    etaNew = etaCoeff(tauRiseNew,tauDecayNew,O.dt/O.superResolution);
    C=filter(1,[etaNew -etaNew*(gammaNew+deltaNew),etaNew*deltaNew],N);  % Update C.
    if O.est_b;
        bNew = mean(Fluorescence-C(O.superResolution:O.superResolution:end));
    else
        bNew = b;
    end;
    C = C + bNew;    
    
    if O.superResolution>1
        cost = BSD_costSR(tauRiseNew,tauDecayNew,Fluorescence,N,bNew,lambda,O.dt,O.superResolution);
    else
        cost =  BSD_cost(tauRiseNew,tauDecayNew,Fluorescence,N,bNew,lambda,O.dt);
    end;

end



%% Refine the offset, and noise amplitude of the signal.


if O.est_sigma ==1 % derivative w.r.t sigma
    MSE=  mean( (Fluorescence-C(O.superResolution:O.superResolution:end)).^2 );
    sigmaNew = sqrt(MSE);    
else
    sigmaNew = sigma;
end;

%% Recompute sparsity prior parameters, Time of tolerance, Threshold, size of signal.
if O.est_a==1
    if O.superResolution>1;
        bias = mean( lambda'./normKernelCoeffSR(tauRise,tauDecay,O.dt,O.superResolution).^2 );
    else
        bias = lambda./normKernelCoeff(tauRise,tauDecay,O.dt)^2;
    end;
    if max(N)>threshold;
        aNew=mean(N(N>threshold)) + bias; % The + takes into account the bias due to L1 regularization, and also acts as a default minimal value of z1/sqrt(normKernelCoeff) for the signal to noise ratio.
    else
        aNew = bias;
    end;
else
    aNew=a;
end;


if O.est_threshold==1
     thresholdNew=thresholdCoeff(tauRiseNew,tauDecayNew,O.dt,sigmaNew,O.z3) ; % Recompute threshold
%    thresholdNew = thresholdCoeff2(tauRiseNew,tauDecayNew,O.dt,sigmaNew,aNew,O.z1,O.z2,O.z3,O.u);
else
    thresholdNew=threshold;
end;


if O.est_lambda==1
    if O.superResolution >1;
%         lambdaNew=sparsityCoeffSR(tauRiseNew,tauDecayNew,O.dt,sigmaNew,O.z1,O.superResolution); % Recompute lambda
        lambdaNew=sparsityCoeffSR2(tauRiseNew,tauDecayNew,O.dt,sigmaNew,aNew,O.z1,O.z2,O.superResolution); % Recompute lambda
    else
%         lambdaNew=sparsityCoeff(tauRiseNew,tauDecayNew,O.dt,sigmaNew,O.z1); % Recompute lambda
        lambdaNew=sparsityCoeff2(tauRiseNew,tauDecayNew,O.dt,sigmaNew,aNew,O.z1,O.z2); % Recompute lambda
    end;
else
    lambdaNew=lambda;
end;


muNew = sparsityCoeffSR(tauRiseNew,tauDecayNew,O.dt,1,1,O.superResolution).^2; % The conditioner for the deconvolution.




end
    
