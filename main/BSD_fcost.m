function cost = BSD_fcost(tauRise,tauDecay,b, sF, sF2,sN,BxN, FxN, NxN,N,mu,dt,T,nmax)

lambdad = exp(-dt/tauDecay);
lambdar = exp(-dt/tauRise);
tmp = 1/((tauDecay/tauRise)^(-(tauRise/(tauDecay-tauRise))) - (tauDecay/tauRise)^(-tauDecay/(tauDecay-tauRise)) );

K = tmp * (lambdad.^[1:nmax+1] - lambdar.^[1:nmax+1])';

K2 = zeros(2*nmax+1,1);
K2(nmax+1:end) = tmp^2 * (lambdad.^[0:nmax]/(1-lambdad^2) - (lambdad.^[0:nmax] + lambdar.^[0:nmax])/(1-lambdad*lambdar) + lambdar.^[0:nmax]/(1-lambdar^2) );
K2(1:nmax) = K2(end:-1:nmax+2);

A = tmp * sum(N(end-nmax:end) .* lambdad.^([nmax:-1:0]') );
B = tmp * sum(N(end-nmax:end) .* lambdar.^([nmax:-1:0]') );


cost = 0.5 * (sF2 + b^2 * T + (NxN' * K2) -2 * ((FxN-b*BxN)' *K) - 2* sF * b  )...
    + mu* sqrt(K2(nmax+1)) * sN - 0.5 * [A;B]' * boundaryCostMatrixCoeff(tauRise,tauDecay,dt) * [A;B];

function matrix = boundaryCostMatrixCoeff(tauRise,tauDecay,dt) 
    matrix = [ [ exp(-dt/tauDecay)^4/(1-exp(-dt/tauDecay)^2), - (exp(-dt/tauDecay) * exp(-dt/tauRise))^2/(1-exp(-dt/tauDecay)*exp(-dt/tauRise))];...
    [- (exp(-dt/tauDecay) * exp(-dt/tauRise))^2/(1-exp(-dt/tauDecay)*exp(-dt/tauRise)), exp(-dt/tauRise)^4/(1-exp(-dt/tauRise)^2)]];
end
end
