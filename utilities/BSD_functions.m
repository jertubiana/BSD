%% A file containing all numerical functions used in the algorithms.

maxiKernelCoeff= @(tauRise,tauDecay) (tauDecay/tauRise)^(-(tauRise/(tauDecay-tauRise))) - (tauDecay/tauRise)^(-tauDecay/(tauDecay-tauRise));

Kernel = @(tauRise,tauDecay,t,dt) max( 1/maxiKernelCoeff(tauRise,tauDecay) * (exp(-(t*dt)/tauDecay) - exp(-(t*dt)/tauRise)), 0);


deltaCoeff = @(tauRise,tauDecay,dt) exp(-dt/tauRise - dt/tauDecay);

gammaCoeff = @(tauRise,tauDecay,dt) exp(-dt/tauRise)+exp(-dt/tauDecay)-exp(-dt/tauRise - dt/tauDecay) ; 

etaCoeff = @(tauRise,tauDecay,dt) (tauRise/tauDecay)^(tauDecay/(tauDecay - tauRise))*( tauDecay/tauRise - 1)./(exp(-dt/tauDecay)-exp(-dt/tauRise));


normKernelCoeff1 = @(tauRise,tauDecay,dt) 1/maxiKernelCoeff(tauRise,tauDecay) * (1/(1-exp(-dt/tauDecay)) - 1/(1-exp(-dt/tauRise))); % sum(K_i)
normKernelCoeff = @(tauRise,tauDecay,dt) 1/maxiKernelCoeff(tauRise,tauDecay) * sqrt(1/(1-exp(-2*dt/tauDecay)) - 2/(1-exp(-dt * (1/tauRise + 1/tauDecay))) + 1/(1-exp(-2*dt/tauRise))); % sqrt(sum(K_i^2))


normKernelCoeffSR1 = @(tauRise,tauDecay,dt,superResolution) 1/maxiKernelCoeff(tauRise,tauDecay) * (...
    exp(-dt/tauDecay).^([superResolution:-1:1]/superResolution)/(1-exp(-dt/tauDecay)) - ...
    exp(-dt/tauRise).^([superResolution:-1:1]/superResolution)/(1-exp(-dt/tauRise))  )';


normKernelCoeffSR = @(tauRise,tauDecay,dt,superResolution) 1/maxiKernelCoeff(tauRise,tauDecay) * sqrt( ...
    exp(-dt/tauDecay).^(2*[superResolution:-1:1]/superResolution)/(1-exp(-dt/tauDecay)^2) -...
     2* (exp(-dt/tauRise)*exp(-dt/tauDecay)).^([superResolution:-1:1]/superResolution)/(1-exp(-dt/tauDecay)*exp(-dt/tauRise)) + ...
    exp(-dt/tauRise).^(2*[superResolution:-1:1]/superResolution)/(1-exp(-dt/tauRise)^2) )';


overlapKernel = @(tauRise,tauDecay,dt,l) 1/maxiKernelCoeff(tauRise,tauDecay)^2 * (...
    exp(-dt/tauDecay).^(2+l)/(1-exp(-dt/tauDecay)^2) -...
    (exp(-dt/tauRise).^l + exp(-dt/tauDecay).^l)*(exp(-dt/tauRise) * exp(-dt/tauDecay))/(1-exp(-dt/tauRise)*exp(-dt/tauDecay)) + ...
    exp(-dt/tauRise).^(2+l)/(1-exp(-dt/tauRise)^2)  );

overlapKernelSR = @(tauRise,tauDecay,dt,l,r1,r2,superResolution) 1/maxiKernelCoeff(tauRise,tauDecay)^2 * (...
    exp(-dt/tauDecay).^(l + (superResolution+1-r1)/superResolution+((superResolution+1-r2)/superResolution))/(1-exp(-dt/tauDecay)^2) - ...
    (exp(-dt/tauDecay).^(l + (superResolution+1-r1)/superResolution) .* exp(-dt/tauRise).^((superResolution+1-r2)/superResolution) + exp(-dt/tauDecay).^((superResolution+1-r2)/superResolution).*exp(-dt/tauRise).^(l + (superResolution+1-r1)/superResolution))/(1-exp(-dt/tauRise)*exp(-dt/tauDecay)) + ...
    exp(-dt/tauRise).^((l + (superResolution+1-r1)/superResolution)+((superResolution+1-r2)/superResolution))/(1-exp(-dt/tauRise)^2) );


boundaryCostMatrixCoeff = @(tauRise,tauDecay,dt) [ [ exp(-dt/tauDecay)^4/(1-exp(-dt/tauDecay)^2), - (exp(-dt/tauDecay) * exp(-dt/tauRise))^2/(1-exp(-dt/tauDecay)*exp(-dt/tauRise))];...
    [- (exp(-dt/tauDecay) * exp(-dt/tauRise))^2/(1-exp(-dt/tauDecay)*exp(-dt/tauRise)), exp(-dt/tauRise)^4/(1-exp(-dt/tauRise)^2)]];


boundaryCostMatrixCoeffSR = @(tauRise,tauDecay,dt,superResolution) [ [ exp(-dt/superResolution/tauDecay)^4/(1-exp(-dt/superResolution/tauDecay)^2), - (exp(-dt/superResolution/tauDecay) * exp(-dt/superResolution/tauRise))^2/(1-exp(-dt/superResolution/tauDecay)*exp(-dt/superResolution/tauRise))];...
    [- (exp(-dt/superResolution/tauDecay) * exp(-dt/superResolution/tauRise))^2/(1-exp(-dt/superResolution/tauDecay)*exp(-dt/superResolution/tauRise)), exp(-dt/superResolution/tauRise)^4/(1-exp(-dt/superResolution/tauRise)^2)]];



sparsityCoeff=@(tauRise,tauDecay,dt,sigma,z1) sigma*normKernelCoeff(tauRise,tauDecay,dt)*z1;
sparsityCoeff2=@(tauRise,tauDecay,dt,sigma,a,z1,z2) min(sigma,normKernelCoeff(tauRise,tauDecay,dt)*a/(z1+z2))*normKernelCoeff(tauRise,tauDecay,dt)*z1;

sparsityCoeffSR = @(tauRise,tauDecay,dt,sigma,z1,superResolution) sigma*normKernelCoeffSR(tauRise,tauDecay,dt,superResolution)*z1;
sparsityCoeffSR2 = @(tauRise,tauDecay,dt,sigma,a,z1,z2,superResolution) BSD_sparsityCoeffSR2(tauRise,tauDecay,dt,sigma,a,z1,z2,superResolution);


% sparsityCoeff_OOPSI = @(sigma,a,nu,dt) sigma^2/(a*nu*dt); 

thresholdCoeff=@(tauRise,tauDecay,dt,sigma,z2) sigma/normKernelCoeff(tauRise,tauDecay,dt)*z2;

thresholdCoeff2=@(tauRise,tauDecay,dt,sigma,a,z1,z2,z3,u) min(...
    u*(a-sparsityCoeff2(tauRise,tauDecay,dt,sigma,a,z1,z2)/normKernelCoeff(tauRise,tauDecay,dt)^2 ),...
    z3 * sigma/normKernelCoeff(tauRise,tauDecay,dt));




