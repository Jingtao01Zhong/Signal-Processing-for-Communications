function [FLOW1, FLOW2] = rw_cfoFilters(FESR,NDEC)
% Get the filters. First one operates at the symbol rate, the second one at
% down-sampled rate. The second filter
    % NDEC = 6;
    % FESR = 1625e3/6 * NDEC;
    rp = 0.1;           % Passband ripple in dB 
    rs = 55;          % Stopband ripple in dB
    f1 = [10e3 270e3];
    dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
    [n,fo,ao,w] = firpmord(f1, [1 0], dev, FESR);
    FLOW1 = firpm(n,fo,ao,w);

    rp = 0.1;           % Passband ripple in dB 
    rs = 40;          % Stopband ripple in dB
    f2 = [3e3 10e3];
    dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
    [n,fo,ao,w] = firpmord(f2, [1 0], dev, FESR/NDEC);
    FLOW2 = firpm(n,fo,ao,w);
end