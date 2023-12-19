function rdsBit = rw_dbpsk_decode(rdsComp,diffMem)
    % Initialize
    rdsBit = zeros(length(rdsComp),1);
    input_bits = zeros(length(rdsComp),1);

    for i = 1 : length(rdsComp)
        if(real(rdsComp(i)) >= 0)
            input_bits(i) = 1;
        end
    end

    % Decode the first bit
    if(diffMem == -1)
        ini_bit = 0;
    else 
        ini_bit = 1;
    end
    if(ini_bit ~= input_bits(1))
        rdsBit(1) = 1;
    end
    % Decode the other bits
    for i = 2 : length(rdsComp)
        if(input_bits(i) ~= input_bits(i-1))
            rdsBit(i) = 1;
        end
    end
    rdsBit = logical(rdsBit);
end