function fltPol = rw_polyint(FPOLY,dec2Pol)
    [num_row, num_col] = size(FPOLY); % number of row is the downsampling factor
    y = zeros(num_row,length(dec2Pol)*num_row);
    for k = 1 : num_row
        y(k,:) = upsample(filter(FPOLY(k,:),1,dec2Pol),num_row,k-1);
    end
    fltPol = sum(y).';
end