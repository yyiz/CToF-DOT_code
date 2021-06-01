function alignedTpsf = alignTpsf_rising_edge(tpsf, refTpsf, inds)
    x = transpose(1:length(tpsf));
    corr = zeros(length(inds), 1);
    
    [~, maxI_ref] = max(refTpsf);
    cropEnd = maxI_ref - 20; % Start cropping at X amount from peak of reference transient
    
    ref_crop = refTpsf(1:cropEnd);
    
    for i = 1:length(inds)
        shift_x = x + inds(i);
        interpMeas = interp1(x, tpsf, shift_x);
        interpMeas(isnan(interpMeas)) = 0;
        
        interp_crop = interpMeas(1:cropEnd);
        
        corr(i) = sum(abs(ref_crop(:)-interp_crop));

    end
    shft = inds(min(corr) == corr);
    if length(shft) > 1
        shft = shft(1);
    end
    shift_x = x + shft;
    alignedTpsf = interp1(x, tpsf, shift_x);
    alignedTpsf(isnan(alignedTpsf)) = 0;
end