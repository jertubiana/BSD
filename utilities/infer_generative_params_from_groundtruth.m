function [Ninf,C,a,b,sigma,tauRise,tauDecay] = infer_generative_params_from_groundtruth(F,N,dt, tauBaseline, tauRise0, tauDecay0,offset)
BSD_functions;

if size(F,1) ==1
    F = F';
end;

if size(N,1) ==1
    N = N';
end;
N = double(N);
F = double(F);


Fluorescence = normalize_remove_baseline(F, 0.15, tauBaseline,dt)


offset = round(offset);
if offset>0;
    N = N(1+offset:end);
    Fluorescence = Fluorescence(1:end-offset);
else
    Fluorescence = Fluorescence(1-offset:end);
    N = N(1:end+offset);
end;

P = struct;
P.tauRise = tauRise0;
P.tauDecay = tauDecay0;
P.b = 0;
O = struct;
O.dt = dt;
O.adaptive = 1;
O.iterations = 200;
O.nNeurons = 1;
O.Time = length(Fluorescence);
[Ninf,C,Palg,~,~]=BSD_ground_truth(N, Fluorescence, O , P);

a = mean(Ninf(N>0));
sigma = std( Fluorescence-C);
b = Palg.b;
tauRise = Palg.tauRise;
tauDecay = Palg.tauDecay;
end


