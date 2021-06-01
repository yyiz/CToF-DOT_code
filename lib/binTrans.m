function binned_trans = binTrans(trans, bf, dim)
    % doesn't include last elements if not evenly divided
    binned_len = floor(size(trans, dim)/bf);
    if dim == 2
        binned_trans = zeros(size(trans,1), binned_len);            
        for ind = 1:binned_len
            binned_trans(:,ind) = sum(trans(:, (bf*(ind-1))+1 : (bf*ind) ), dim);
        end
    elseif dim == 3
        binned_trans = zeros(size(trans,1), size(trans,2), binned_len);
        for ind = 1:binned_len
            binned_trans(:,:,ind) = sum(trans(:, :, (bf*(ind-1))+1 : (bf*ind) ), dim);
        end
    else
        error('Dim not 2 or 3');
    end
    
end