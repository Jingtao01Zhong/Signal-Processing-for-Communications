function fmSig = rw_fmrx(rxSig) % Samples to read at one round is 4096*4
    sectionLength = length(rxSig);
    delay_factor = 1;

    path_delay = rxSig(1:sectionLength-delay_factor);
    path_conjunction = conj(rxSig(1+delay_factor:sectionLength));
    phase_detect = path_conjunction .* path_delay;
    fmSig = angle(phase_detect);
    fmSig = fmSig - mean(fmSig);
end