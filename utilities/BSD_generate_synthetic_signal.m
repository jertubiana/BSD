function [N,C,F] = BSD_generate_synthetic_signal(P,O,varargin)

BSD_functions;

P.delta=deltaCoeff(P.tauRise,P.tauDecay,O.dt);
P.gamma=gammaCoeff(P.tauRise,P.tauDecay,O.dt);
P.eta=etaCoeff(P.tauRise,P.tauDecay,O.dt);

if ~isempty(varargin)
    N = varargin{1};
    if size(N,1) ~= O.Time;
        N = N';
    end;
else
    if isfield(O,'isolated');
        tau = round(20* (P.tauRise+P.tauDecay)/O.dt);
        N = zeros(O.Time,O.nNeurons);
        indexes = 1:tau:O.Time;
        if isfield(O,'superResolution');
            n =O.superResolution;
        else
            n = 20;
        end;
        indexes = indexes + randi(n,1,length(indexes));
        indexes = min(indexes,O.Time);
        N(indexes,:) = 1;
    else
        N = poissrnd(P.nu*O.dt*ones(O.Time,O.nNeurons)); % simulate spike train
    end;
end;
    

noise = randn(O.Time,O.nNeurons); % noise

C = zeros(size(N));
F = zeros(size(N));
for i = 1:O.nNeurons;
    C(:,i) = P.a* filter(1,[1 -P.gamma-P.delta,P.delta],N(:,i))/P.eta + P.b;
    F(:,i) = C(:,i) + P.sigma * noise(:,i);
end;
end
