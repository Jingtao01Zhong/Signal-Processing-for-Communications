function window = windowSelect(windowType, windowLength)

switch windowType

    case 'barthannwin'
        window = barthannwin(windowLength); % Modified Bartlett-Hann window
    case 'bartlett'
        window = bartlett(windowLength);
    case 'blackman'
        window = blackman(windowLength); 
    case 'blackmanharris'
        window = blackmanharris(windowLength);
    case 'bohmanwin'
        window = bohmanwin(windowLength); 
    case 'chebwin'
        window = chebwin(windowLength); 
    case 'enbw'
        window = enbw(windowLength); 
    case 'flattopwin'
        window = flattopwin(windowLength); 
    case 'gausswin'
        window = gausswin(windowLength); 
    case 'hamming'
        window = hamming(windowLength);
    case 'hann'
        window = hann(windowLength); 
    case 'kaiser'
        window = kaiser(windowLength); 
    case 'nuttallwin'
        window = nuttallwin(windowLength); 
    case 'parzenwin'
        window = parzenwin(windowLength);
    case 'rectwin'
        window = rectwin(windowLength); 
    case 'taylorwin'
        window = taylorwin(windowLength); 
    case 'triang'
        window = triang(windowLength);
    case 'tukeywin'
        window = tukeywin(windowLength);
    otherwise
        error('please enter the correct window name');
end

end
