function [burst,flg] = rw_getBurst(rxFlt2)
    nSym_FCCH = 142;
    burst = zeros(1,3*nSym_FCCH);

    for k = 1:3 % get 3 burst signals
        rxFlt2_cut = rxFlt2(1 + (k-1)*20*nSym_FCCH : 40*nSym_FCCH + (k-1)*20*nSym_FCCH);
        for i = 1 : 19*nSym_FCCH % get at least 2 burst signals
            % Calculate the energy
            E(i) = sum(abs(rxFlt2_cut(i:i+nSym_FCCH-1)).^2);
        end
        
        [~,i_max] = max(E);
        burst((k-1)*nSym_FCCH+1:k*nSym_FCCH) = rxFlt2_cut(i_max :i_max+nSym_FCCH-1);
        
        % find the next burst signal
        for i = 1:10*nSym_FCCH
            E_nxt(i) = sum(abs(rxFlt2_cut(i+i_max+5*nSym_FCCH/2:i+nSym_FCCH-1+i_max+5*nSym_FCCH/2)).^2); 
        end
        [~,i_submax] = max(E_nxt);
        i_submax = i_submax + i_max+5*nSym_FCCH/2;
    
        % calculate the difference of the burst indices
        d = i_submax - i_max;

        if(d < 8*nSym_FCCH || d > 12*nSym_FCCH)
            detect(k) = 1;
        else 
            detect(k) = 0;
        end
    end
    flg = ismember(1, detect);
end