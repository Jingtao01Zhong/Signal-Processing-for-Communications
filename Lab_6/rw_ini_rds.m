function [FRDS,MANMF] = rw_ini_rds(FESR,L)
    % FESR = 228e3， NDEC = 12
    rp = 1;           % Passband ripple in dB 
    rs = 40;          % Stopband ripple in dB
    f1 = [1e3 10e3];
    dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
    [n,fo,ao,w] = firpmord(f1, [1.2 0], dev, FESR);
    FRDS = firpm(n,fo,ao,w);

    beta = 1; % roll off factor
    span = 3;    % truncated to SPAN symbols and each symbol is represented by SPS samples
    sps = L/2;     % samples per symbol
    h = rcosdesign(beta, span, sps);

    delay = L/4;
    delayed_h = [zeros(1, delay), h(1:end)];
    h = [h, zeros(1, delay)];
    MANMF = h - delayed_h;
end