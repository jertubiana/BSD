function DeltaF = normalize_remove_baseline(F, q, tauBaseline,dt)
DeltaF = F;
DeltaF = (DeltaF-median(DeltaF))/std(DeltaF);
window = min(round(tauBaseline/dt),length(F));
baseline0 = running_percentile(DeltaF,window, q * 100)';
DeltaF = DeltaF - baseline0;