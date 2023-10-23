% rw_window is a function that returns the low-pass filter coefficients.
function coefficients = rw_window(fs) % fs is the sampling frequency

    %hamming window
    wp = 20e3 * pi / (fs/2); % passband
    ws = 34e3 * pi / (fs/2); % stopband
    wc = (wp+ws)/2;
    M = round(8*pi/(ws-wp));
    n = 0 : M;
    h=wc/pi*sinc(wc/pi*(n-M/2)/2);
    func_window = h.' .* hamming(M+1);
    coefficients = func_window;
    % HH = freqz(func_window,1);
    % ww = linspace(0,fs/2,length(HH));
    % plot(ww,20*log10(abs(HH))), ylabel Magnitude, xlabel Frequency
    % title('Magnitude response')
end