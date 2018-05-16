function xshuff = phaseperm(x)
xfft = fft(x);
ph = imag(xfft);
xshuff = abs(ifft(real(xfft).*exp(sqrt(-1)*ph(randperm(length(ph))))));
