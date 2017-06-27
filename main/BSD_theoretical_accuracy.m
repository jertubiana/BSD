function [RF,DeltaT,lag,width,FPR,FNR] = BSD_theoretical_accuracy(P,O)



if ~isfield(O,'z1'); O.z1 = norminv(1-1e-2,0,1); end;
if ~isfield(O,'z2'); O.z2 = norminv(1-1e-2,0,1); end;
if ~isfield(O,'delta_max'); O.delta_max = 20;  end;
if ~isfield(O,'Nsim'); O.Nsim = 1e4; end;
if ~isfield(O,'superResolution'); O.superResolution = 1; end;
if ~isfield(O,'discretisation'); O.discretisation = 40; end;


n = P.sigma/(P.a*normKernelCoeff(P.tauRise,P.tauDecay,O.dt));

z = min(O.z1,O.z1/(n*(O.z1+O.z2))); % The normalized sparsity parameter.

if O.superResolution >0;
    O.discretisation = ceil(O.discretisation/O.superResolution) * O.superResolution;
    [RF,DeltaT] = probability_mismatch_SR(P.tauRise,P.tauDecay,O.dt,P.sigma,P.a,z,O.superResolution,O.delta_max,O.discretisation,O.Nsim);
else
    [RF,DeltaT] = probability_mismatch(P.tauRise,P.tauDecay,O.dt,n,z,O.delta_max,O.discretisation,O.Nsim);
end;

options = optimoptions(@lsqcurvefit,'Display','off');
bell = @(params,X) params(1)*exp(-(X-params(2)).^2/(2*params(3)^2));
x = lsqcurvefit(bell,[1,0,O.dt],DeltaT,RF,[0,-Inf,0],[],options);
lag = x(2);
width = x(3);

FNR = 1-sum(RF); % False negative rate.

FPR = (1-normcdf(z))/O.dt; % # Fake non-zero spikes/s.

RF = RF/(DeltaT(2)-DeltaT(1)); % Make it a probability density.

function [proba,DeltaT] = probability_mismatch(tauRise,tauDecay,dt,n,z,delta_max,discretisation,Nsim)
proba = zeros(2*discretisation*delta_max+1,1);
norm = normKernelCoeff(tauRise,tauDecay,dt);
overlaps = overlapKernel(tauRise,tauDecay,dt,[0:2*delta_max]);
KtK = toeplitz(overlaps)/norm^2;
generative = KtK^0.5;
epsilon =  randn(Nsim,2*delta_max+1) * generative;

for r = 1:discretisation;
    overlap_signal_left = overlapKernelSR(tauRise,tauDecay,dt,[1:delta_max],1,r,discretisation);
    overlap_signal_right = overlapKernelSR(tauRise,tauDecay,dt,[0:delta_max],r,1,discretisation);
    overlap_signal = [overlap_signal_left(end:-1:1),overlap_signal_right]/norm^2;
    display(overlap_signal*norm^2);
    lik = -0.5*  max(bsxfun(@plus,n*epsilon, overlap_signal -n*z),0).^2;
    [ val,pos_min] = min( lik, [], 2 );
    undetected = (val == 0);

    for i = 1:2*delta_max+1
        bool = (pos_min==i) & ~undetected;
        if (sum(bool)>0);
            display(sum(bool))
%             display([i,r]);
            display([i * discretisation + (1-r)]);
        end;
        I = [max((i-1) * discretisation + (2-r),1): min(i * discretisation + (1-r),2*delta_max*discretisation+1)];
        proba(I) = proba(I) + sum(bool)/Nsim/discretisation^2;
    end
end;
DeltaT = [-delta_max * discretisation: delta_max * discretisation]' * dt/discretisation;
end


function [proba,DeltaT] = probability_mismatch_SR(tauRise,tauDecay,dt,sigma,a,z,superResolution,delta_max,discretisation,Nsim)
proba = zeros(2*discretisation*delta_max+1,1);

norms = normKernelCoeffSR(tauRise,tauDecay,dt,superResolution);

% % % Generate the noise
nmax = 2*delta_max;
[ls,r1,r2] = ndgrid( [0:nmax], [1:superResolution], [1:superResolution] );
tmp = zeros(2*nmax+1,superResolution,superResolution);
tmp(nmax+1:end,:,:) = overlapKernelSR(tauRise,tauDecay,dt,ls,r1,r2,superResolution);


for k = 1:nmax;
    tmp(k,:,:) = permute(tmp(2*(nmax+1)-k,:,:),[1,3,2]);
end;

autocov = zeros((delta_max+1) * superResolution);
for Delta = [-delta_max:delta_max];
    for DeltaP = [-delta_max:delta_max];
        index = (DeltaP - Delta) + nmax+1;
        autocov( (delta_max+Delta) * superResolution + [1:superResolution],...
            (DeltaP+delta_max)* superResolution+[1:superResolution])  = ...
            tmp(index,:,:);
    end;
end;

[V,Lam] = eig(autocov);
index = find(diag(Lam)>1e-6);
rank = length(index);
generative = sqrt(Lam(index,index)) * V(:,index)';
epsilon =  reshape(randn(Nsim,rank) * generative,[Nsim,2*delta_max+1,superResolution]);
% % % Compute the overlaps.

step = discretisation/superResolution;

for r = 1:discretisation;
    [Deltas, rs] = ndgrid([1:delta_max], [1:step:discretisation]); % Reconstruction indexes.
    overlap_signal_left = overlapKernelSR(tauRise,tauDecay,dt,Deltas,rs,r,discretisation);
    [Deltas, rs] = ndgrid([0:delta_max], [1:step:discretisation]); % Reconstruction indexes.
    
    overlap_signal_right = overlapKernelSR(tauRise,tauDecay,dt,Deltas,r,rs,discretisation);
    overlap_signal = [overlap_signal_left(end:-1:1,:);overlap_signal_right];
    lik = -0.5 * bsxfun(@rdivide,...
        max(...
        bsxfun(@plus, sigma/a * epsilon, ...
        reshape(...
        bsxfun(@minus,overlap_signal,z * sigma/a * norms')...
        ,[1,2*delta_max+1,superResolution])...
        ) ...
        ,0).^2, ...
        reshape(norms.^2,[1,1,superResolution])...
        ); 
    
    [ val,pos_min] = min( reshape( permute( lik,[1,3,2]),[Nsim,superResolution * (2*delta_max+1)]), [], 2 ); % Care, fortran ordering of array in matlab...
    undetected = (val == 0);
    for l = 1: (2*delta_max+1)*superResolution;
        bool = (pos_min==l) & ~undetected;
        I =  [max( (l-1) * step + (2-r),1): min( l * step + (1-r),2*delta_max*discretisation+1)];
        proba(I) = proba(I) + sum(bool)/Nsim/(discretisation*step);
%         if sum(bool)>0;
%             display(sum(bool))
%             display(l);
%             display(r);
%             plot(proba);            
%         end;        
    end;
end    
DeltaT = [-delta_max * discretisation: delta_max * discretisation]' * dt/discretisation;    
end


    

function M = maxiKernelCoeff(tauRise,tauDecay)
    M = (tauDecay/tauRise)^(-(tauRise/(tauDecay-tauRise))) - (tauDecay/tauRise)^(-tauDecay/(tauDecay-tauRise));
end


function N = normKernelCoeff(tauRise,tauDecay,dt)
    N = sqrt((exp(-dt/tauDecay)-exp(-dt/tauRise))^2*(1+exp(-dt/tauRise-dt/tauDecay))...
    /((1-exp(-2*dt/tauRise))*(1-exp(-2*dt/tauDecay))*(1-exp(-dt/tauRise-dt/tauDecay))))*(tauDecay/((tauDecay-tauRise)*(tauRise/tauDecay)^(tauRise/(tauDecay-tauRise))));
end

function N = normKernelCoeffSR(tauRise,tauDecay,dt,superResolution)
    N = 1/maxiKernelCoeff(tauRise,tauDecay) * sqrt( ...
    exp(-dt/tauDecay).^(2*[superResolution:-1:1]/superResolution)/(1-exp(-dt/tauDecay)^2) -...
     2* (exp(-dt/tauRise)*exp(-dt/tauDecay)).^([superResolution:-1:1]/superResolution)/(1-exp(-dt/tauDecay)*exp(-dt/tauRise)) + ...
    exp(-dt/tauRise).^(2*[superResolution:-1:1]/superResolution)/(1-exp(-dt/tauRise)^2) )';
end


function overlap = overlapKernel(tauRise,tauDecay,dt,l)
    overlap = 1/maxiKernelCoeff(tauRise,tauDecay)^2 * (...
    exp(-dt/tauDecay).^(2+l)/(1-exp(-dt/tauDecay)^2) -...
    (exp(-dt/tauRise).^l + exp(-dt/tauDecay).^l)*(exp(-dt/tauRise) * exp(-dt/tauDecay))/(1-exp(-dt/tauRise)*exp(-dt/tauDecay)) + ...
    exp(-dt/tauRise).^(2+l)/(1-exp(-dt/tauRise)^2)  );
end

function overlapSR = overlapKernelSR(tauRise,tauDecay,dt,l,r1,r2,superResolution)

overlapSR = 1/maxiKernelCoeff(tauRise,tauDecay)^2 * (...
    exp(-dt/tauDecay).^(l + (superResolution+1-r1)/superResolution+((superResolution+1-r2)/superResolution))/(1-exp(-dt/tauDecay)^2) - ...
    (exp(-dt/tauDecay).^(l + (superResolution+1-r1)/superResolution) .* exp(-dt/tauRise).^((superResolution+1-r2)/superResolution) + exp(-dt/tauDecay).^((superResolution+1-r2)/superResolution).*exp(-dt/tauRise).^(l + (superResolution+1-r1)/superResolution))/(1-exp(-dt/tauRise)*exp(-dt/tauDecay)) + ...
    exp(-dt/tauRise).^((l + (superResolution+1-r1)/superResolution)+((superResolution+1-r2)/superResolution))/(1-exp(-dt/tauRise)^2) );
end

end
