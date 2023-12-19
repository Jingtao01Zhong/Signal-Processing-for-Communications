function [lpSig1,lpSig2] = rw_bank(fmSig, FLOW, FESR)
    % Calculate FFT of fmSig and FLOW
    N = length(fmSig) + length(FLOW) - 1;
    n = ceil(N/8); % NDEC = 8
    N = n*8;
    fmSig_FFT = fft(fmSig,N);
    FLOW_FFT = fft(FLOW.',N);

    % Make the filterbank: the FFT of FLOW + its circular shift in frequency domain
    shift = ceil(38e3 / FESR * N);
    Filter_Bank = FLOW_FFT + circshift(FLOW_FFT,shift);

    % Multiply the FFT of fmSig with the filter bank
    y = fmSig_FFT .*  Filter_Bank;
    % Extract the two channels. This makes decimation in frequency domain.
    mono_part = floor([1 : 15e3/FESR*N , 225e3/FESR*N+1 : N]); 
    stereo_part = floor([30e3/FESR*N+1 : 53e3/FESR*N , 23e3/FESR*N : 30e3/FESR*N]);

    % Take ifft() of the L+R and L-R channels and make the output real() to take care
    % of the numerical errors
    sig_mono = ifft(y(mono_part));
    sig_stereo = ifft(y(stereo_part));
    lpSig1 = real(sig_mono);
    lpSig2 = real(sig_stereo);
    % Return the two FM channels in time domain. The sampling rate in each
    % channel is 30 KHz. Column vectors.
    
end