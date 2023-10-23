function [FPASS, FPROTO, FLOW, FPROTOI, FINT, FPER, NDEC, FESR, FPOLY]  = rw_ini_filt()
    FESR = 300e3; % define the samping rate
    NDEC = 15; % define the down-samping factor

    % design FPASS (direct implementation and the prototype filter)
    rp = 0.5;           % Passband ripple in dB 
    rs = 65;          % Stopband ripple in dB
    dev = [10^(-rs/20) (10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
    f = [15e3 18e3 20e3 23e3]; % stopband and passband frequency
    [n,fo,ao,w] = firpmord(f, [0 1 0], dev, FESR);
    FPASS = firpm(n,fo,ao,w);

    % design FPROTO for the polyphase filter
    f = [1e3*19/4 4e3*19/4]; % passband and stopband frequency
    dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
    [n,fo,ao,w] = firpmord(f, [1 0], dev, FESR);
    FPROTO = firpm(n,fo,ao,w);

    % generate the matrix whose rows contain the polyphase components
    FPOLY = zeros(NDEC,ceil(length(FPROTO)/NDEC));
    for i = 1 : NDEC
        x = FPROTO(i:NDEC:end);
        if(length(x) < ceil(length(FPROTO))/NDEC)
            FPOLY(i,:) = [x 0];
        else 
            FPOLY(i,:) = x;
        end
    end

    % FLOW: Low-pass filter running at down-sampled frequency to extract the carrier
    f = [1e3 3e3];
    [n,fo,ao,w] = firpmord(f, [1 0], dev, FESR/NDEC);
    FLOW = firpm(n,fo,ao,w);

    % FPROTOI: prototype filter
    f = [20e3 60e3];
    rp = 0.5;           % Passband ripple in dB 
    rs = 55;          % Stopband ripple in dB
    dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
    [n,fo,ao,w] = firpmord(f, [1 0], dev, FESR);
    FPROTOI = firpm(n,fo,ao,w);

    % FPER: sparse version with multiple compressed copies of FPROTOI's frequency response
    FPER = zeros(1,length(FPROTOI)*NDEC);
    for i = 1: length(FPROTOI)
        FPER((i-1)*NDEC + 1) = FPROTOI(i);
    end

    % FINT: Interpolating filter
    rp = 0.5;           % Passband ripple in dB 
    rs = 55;          % Stopband ripple in dB
    dev = [10^(-rs/20) (10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
    f = [1e3 18e3 20e3 39e3]; % stopband and passband frequency
    [n,fo,ao,w] = firpmord(f, [0 1 0], dev, FESR);
    FINT = firpm(n,fo,ao,w);

end