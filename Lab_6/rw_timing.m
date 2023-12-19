function rdsSym = rw_timing(rdsMatched,L)
    % Energy Detector
    for i = 1:L
        E(i) = sum(abs(rdsMatched(i:L:end)).^2);
    end
    [~,i_max] = max(E);
    rdsSym = rdsMatched(i_max:L:end);
end