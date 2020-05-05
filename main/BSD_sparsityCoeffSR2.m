function lambda = BSD_sparsityCoeffSR2(tauRise,tauDecay,dt,sigma,a,z1,z2,superResolution)


norm = normKernelCoeffSR(tauRise,tauDecay,dt,superResolution);

sigma = max(a/100,sigma); % If superResolution > 1, add L1 penalty to secure optimization convergence.
fac = min(sigma, a * mean(norm)/(z1+z2) );
lambda = z1 * norm *fac;


function norm = normKernelCoeffSR(tauRise,tauDecay,dt,superResolution) 
	M = (tauDecay/tauRise)^(-(tauRise/(tauDecay-tauRise))) - (tauDecay/tauRise)^(-tauDecay/(tauDecay-tauRise));
	norm = 1/M * sqrt( ...
    exp(-dt/tauDecay).^(2*[superResolution:-1:1]/superResolution)/(1-exp(-dt/tauDecay)^2) -...
     2* (exp(-dt/tauRise)*exp(-dt/tauDecay)).^([superResolution:-1:1]/superResolution)/(1-exp(-dt/tauDecay)*exp(-dt/tauRise)) + ...
    exp(-dt/tauRise).^(2*[superResolution:-1:1]/superResolution)/(1-exp(-dt/tauRise)^2) );
end


end

