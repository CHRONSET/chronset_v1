function [goP] = compute_goodness_of_pitch(fdt);
%%
goP = zeros(size(fdt,1),1);
for jt = 1:size(fdt,1)
    
    x = xcorr(log(abs(fdt(jt,:))),'coef');
    
    nfft = 2^nextpow2(length(x));
    y = fft(x,nfft);
    y = real(y).^2+imag(y).^2;
    y = y(1:nfft/2+1);
    
    
    goP(jt) = max(y);
    
end;