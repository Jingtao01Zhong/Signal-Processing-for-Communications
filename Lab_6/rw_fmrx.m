function fmSig = rw_fmrx(rxSig)
    sectionLength = length(rxSig);
    delay_factor = 1;
    
    path_delay = [rxSig(1:sectionLength-delay_factor);0];
    path_conjunction = conj([rxSig(1+delay_factor:sectionLength);0]);
    phase_detect = path_conjunction .* path_delay;
    fmSig = angle(phase_detect);
    fmSig = fmSig - mean(fmSig);
end