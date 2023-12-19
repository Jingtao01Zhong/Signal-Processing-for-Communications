function FLOW = rw_ini_bank()
% Design the protype filter FLOW using least-squares. Using fir2() is fine. 
% The order of the filter must be a multiple of the down-sampling factor and the
% attenuation at 19 KHz at least 40dB
    f = [0 14e3 17.5e3 240e3/2];
    f = f/(240e3/2);
    m = [1 1 0 0];
    FLOW = fir2(8*20,f,m);
end