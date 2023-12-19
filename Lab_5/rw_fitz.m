function fitzoffs = rw_fitz(burst, fsample, MAXLAG)
% length(burst) = 142
% nsample = FESR/NDEC = 2.7803e05
% MAXLAG = 5

    T = 1/fsample;
    % M = 1/2/T/f_max; M <= MAXLAG = 5
    % f_max = 1/2/MAXLAG/T;

    N = length(burst);
    f_est = 0;
    for i = 1:MAXLAG
        M = i;
        for m = 1 : M
            R_sum = 0;
            for n = m+1 : N
                R_sum = R_sum + burst(n) * conj(burst(n-m));
            end
            R(m) = R_sum/(N-m);
            f_est = f_est + angle(R(m)); 
        end
        fitzoffs(i) = f_est/(pi*M*(M+1)*T);
    end
end
