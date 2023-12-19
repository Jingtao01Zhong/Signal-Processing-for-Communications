function rdsComp = rw_offset(rdsSym,NV)
% Viterbi&Viterbi phase offset estimation. The second argument is the
% length of the moving average. The output is the signal where the
% phase offset is compensated
% NV = 24
    M = 2; % BPSK
    sig = rdsSym.^M;
    sig = sig./abs(sig);
    H_avg = ones(1, NV+1) / (NV+1);

    sig_avg_r = filter(H_avg,1,real(sig));
    sig_avg_i = filter(H_avg,1,imag(sig));
    sig_avg = sig_avg_r + 1j*sig_avg_i;

    % sig_avg = filter(H_avg,1,sig);

    sig_angle = angle(sig_avg)/M;
    
    N = 1;
    y = zeros(length(rdsSym),1);
    y(1:end-N+1) = rdsSym(N:end);
    y(end-N+1:end) = rdsSym(end);

    rdsComp = exp(-1j * sig_angle) .* y;
end