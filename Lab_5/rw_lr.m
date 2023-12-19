function lroffs = rw_lr(burst, fsample, MAXLAG)
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
            f_est = f_est + R(m); 
        end
        lroffs(i) = 1/(pi*(M+1)*T) * angle(f_est);
    end
end