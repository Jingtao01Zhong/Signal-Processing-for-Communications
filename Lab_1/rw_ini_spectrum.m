function [nSection, nOverlap, wintype, nFFT] = rw_ini_spectrum()

nSection = 8;
nOverlap = 0.4;
%wintype = 'rectwin'; %kaiser
wintype = "hann";
nFFT = 1024;

end