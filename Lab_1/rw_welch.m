%% Own spectrum analyzer
function perio = rw_welch(rxSig, nSection, nOverlap, wintype, nFFT)

%calculate the length of segement
%nSection = (dataLength-nOverlap*sectionLength)/[sectionLength*(1-nOverlap)]
sectionLength = floor(length(rxSig)/nSection/(1-nOverlap+nOverlap/nSection));

%initialize the result vector
spectrumResult = zeros(1, nFFT);
windowLength = sectionLength;
window = windowSelect(wintype, windowLength);
% calculate the fft and the power 
for i = 1: nSection
    %calculate the start Index and end Index
    startIndex = floor((i - 1) * sectionLength * (1 - nOverlap) + 1);
    endIndex = startIndex + sectionLength - 1;

    %calculate the PSD
    section = rxSig(startIndex:endIndex) .* window;
    fftResult = fftshift(fft(section, nFFT));
    powerSpectrum = 1/sectionLength * abs(fftResult).^2;
    spectrumResult = spectrumResult + powerSpectrum.';
end

spectrumResult = spectrumResult / nSection;
perio = spectrumResult.';