function fm57 = rw_shift_freq(fmSig,FESR,f_demodu)
    % f_demodu = -57e3
    t = 0 : 1/FESR : (length(fmSig)-1)/FESR;
    signal_demodu = exp(1j * 2*pi*f_demodu*t);
    fm57 = signal_demodu.' .* fmSig;
end