function N_resampled = resample_spike(N,f,f_target)

s = size(N);
if s(2) ~=1;
    N = N';
end;
    

ratio = f_target/f;
if ratio<1
    N_resampled = zeros(length(N)*ratio,1);
    for k = 1:length(N);
        id_resampled = 1 + floor( (k-1)*ratio);
        N_resampled(id_resampled ) =  N_resampled(id_resampled) + N(k);
    end;
elseif ratio>1
    N_resampled = zeros(ratio* length(N), 1);
    for k=1:length(N);
        N_resampled(1+(k-1)*ratio:k*ratio) = N(k);%/ratio;
    end;
else
    N_resampled = N;
end;

if s(2) ~=1;
    N_resampled = N_resampled';
end;