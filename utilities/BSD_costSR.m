function cost = BSD_costSR(tauRise,tauDecay,Fluorescence,N,b,lambda,dt,superResolution)

BSD_functions;
delta=deltaCoeff(tauRise,tauDecay,dt/superResolution);
gamma=gammaCoeff(tauRise,tauDecay,dt/superResolution);
eta=etaCoeff(tauRise,tauDecay,dt/superResolution);

C=filter(1,[eta, -eta * (gamma+delta)  ,  eta * delta  ],N);
downC = C(superResolution:superResolution:end);
error=Fluorescence-downC-b;
sN = zeros(superResolution,1);
for i =1:superResolution; sN(i) = sum(N(i:superResolution:end)); end;

cost = 0.5*sum(error.^2)+ sum(lambda' .* sN);
end