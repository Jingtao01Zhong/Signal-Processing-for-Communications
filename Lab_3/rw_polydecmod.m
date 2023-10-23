function decPol = rw_polydecmod(FPOLY,rxSig,freqZone) % freqZone = 1
    [num_row, num_col] = size(FPOLY); % number of row is the downsampling factor
    y = zeros(num_row,length(rxSig)/num_row);
    for k = 1 : num_row
        y(k,:) = filter(FPOLY(k,:),1,rxSig(k:num_row:end)) .* (exp(2*pi*(k-1)*(freqZone)/num_row*1i)); %exp(j*2*pi/NDEC*[0:nSample-1].')
    end
    decPol = sum(y).';
end