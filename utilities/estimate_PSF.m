function PSF = estimate_PSF(N,Ninf,maxlag)
% NxN = xcov(N,N,2*maxlag,'biased');
% NinfxN = xcov(N,Ninf,maxlag,'biased');
NxN = xcorr(N,N,2*maxlag);
NinfxN = xcorr(N,Ninf,maxlag);
PSF = toeplitz(NxN(1+2*maxlag:end))\NinfxN;
PSF = PSF(end:-1:1);
