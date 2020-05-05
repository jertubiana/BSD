function cost = BSD_fcostSR(tauRise,tauDecay,b, sF, sF2,sN,BxN, FxN, NxN, N, mu, dt,T,nmax,superResolution)

lambdad = exp(-dt/tauDecay);
lambdar = exp(-dt/tauRise);
M = (tauDecay/tauRise)^(-(tauRise/(tauDecay-tauRise))) - (tauDecay/tauRise)^(-tauDecay/(tauDecay-tauRise));

[ls, rho] = ndgrid([0:nmax],[superResolution:-1:1]/superResolution);
kernel = max( lambdad.^(ls + rho) - lambdar.^(ls+rho) , 0)/M;



[ls,r1,r2] = ndgrid( [0:nmax], [1:superResolution], [1:superResolution] );
kernel_autocov = zeros(2*nmax+1,superResolution,superResolution);
kernel_autocov(nmax+1:end,:,:) = overlapKernelSR(tauRise,tauDecay,dt,ls,r1,r2,superResolution);

for k = 1:nmax;
    kernel_autocov(k,:,:) = kernel_autocov(2*(nmax+1)-k,:,:);
    kernel_autocov(2*(nmax+1)-k,:,:) = permute(kernel_autocov(2*(nmax+1)-k,:,:),[1,3,2]);
end;


kernel_autocov = kernel_autocov/M^2;
normK2 = sqrt(diag(squeeze(kernel_autocov(nmax+1,:,:))) );

lambdas = mu .* normK2;



error =  0.5 * (sF2 +  sum(kernel_autocov(:) .* NxN(:)) + b.^2 * T)...
    - sum(kernel(:) .* (FxN(:)- b*BxN(:)))  - sF * b;


A = 1/M * sum(N(end-nmax*superResolution:end) .* lambdad.^([nmax*superResolution:-1:0]'/superResolution) );
B = 1/M * sum(N(end-nmax*superResolution:end) .* lambdar.^([nmax*superResolution:-1:0]'/superResolution) );


cost = error + sum(lambdas .* sN) - 0.5 * [A;B]' * boundaryCostMatrixCoeffSR(tauRise,tauDecay,dt,superResolution) * [A;B];

function overlap = overlapKernelSR(tauRise,tauDecay,dt,l,r1,r2,superResolution)
    overlap =  (...
    exp(-dt/tauDecay).^(l + (superResolution+1-r1)/superResolution+((superResolution+1-r2)/superResolution))/(1-exp(-dt/tauDecay)^2) - ...
    (exp(-dt/tauDecay).^(l + (superResolution+1-r1)/superResolution) .* exp(-dt/tauRise).^((superResolution+1-r2)/superResolution) + exp(-dt/tauDecay).^((superResolution+1-r2)/superResolution).*exp(-dt/tauRise).^(l + (superResolution+1-r1)/superResolution))/(1-exp(-dt/tauRise)*exp(-dt/tauDecay)) + ...
    exp(-dt/tauRise).^((l + (superResolution+1-r1)/superResolution)+((superResolution+1-r2)/superResolution))/(1-exp(-dt/tauRise)^2) );
end

function matrix = boundaryCostMatrixCoeffSR(tauRise,tauDecay,dt,superResolution) 
    matrix = [ [ exp(-dt/superResolution/tauDecay)^4/(1-exp(-dt/tauDecay)^2), - (exp(-dt/superResolution/tauDecay) * exp(-dt/superResolution/tauRise))^2/(1-exp(-dt/tauDecay)*exp(-dt/tauRise))];...
    [- (exp(-dt/superResolution/tauDecay) * exp(-dt/superResolution/tauRise))^2/(1-exp(-dt/tauDecay)*exp(-dt/tauRise)), exp(-dt/tauRise)^4/(1-exp(-dt/tauRise)^2)]];
end


end
