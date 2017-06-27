function cost = BSD_cost(tauRise,tauDecay,F,N,b,lambda,dt)

BSD_functions;
delta=deltaCoeff(tauRise,tauDecay,dt);
gamma=gammaCoeff(tauRise,tauDecay,dt);
eta=etaCoeff(tauRise,tauDecay,dt); 

C=filter(1,[eta, -eta * (gamma+delta)  ,  eta * delta  ],N);
error=F-C-b;
cost = 0.5*sum(error.^2)+ lambda *sum(N);


end