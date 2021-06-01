
%% 
mS = m_abs;
mB = m_bkg;

% reduce dimensionality by gating and binning
gate = 3001:5000;
bin_factor = 4;
mS = binTrans(mS(:,gate),bin_factor,2);
mB = binTrans(mB(:,gate),bin_factor,2);
bJac = binTrans(bkgTransients(:,:,gate),bin_factor,3);
sJac = binTrans(absTransients(:,:,gate),bin_factor,3);
refTrans = mB(1,:);

% align jacobian to bkg reference
inds = -25:.1:25; 
fprintf('\nAligning transient: %04d/%04d',0,size(mS,1));
counter = 1;
for r = 1:size(sJac, 1)
    for c = 1:size(sJac, 2)
        fprintf('\b\b\b\b\b\b\b\b\b%04d/%04d',counter,size(sJac,1)*size(sJac,2));
        sJac(r,c,:) = alignTpsf_rising_edge(...
            squeeze(sJac(r,c,:)),...
            refTrans, inds);
        bJac(r,c,:) = alignTpsf_rising_edge(...
            squeeze(bJac(r,c,:)),...
            refTrans, inds);
        counter = counter + 1;
    end
end
fprintf('\nFinished aligning.\n');

% align measurements to bkg reference
fprintf('\nAligning transient: %03d/%03d',0,size(mS,1));
for ind = 1:size(mS, 1)
    fprintf('\b\b\b\b\b\b\b%03d/%03d',ind,size(mS,1));
    mS(ind,:) = alignTpsf_rising_edge(...
        mS(ind,:),...
        refTrans, inds);
    mB(ind,:) = alignTpsf_rising_edge(...
        mB(ind,:),...
        refTrans, inds);
end
fprintf('\nFinished aligning.\n');

%% apply further gating
gate2 = 1:500;
bin_factor2 = 250;
mS_g = binTrans(mS(:,gate2),bin_factor2,2);
mB_g = binTrans(mB(:,gate2),bin_factor2,2);
bJac_g = binTrans(bJac(:,:,gate2),bin_factor2,3);
sJac_g = binTrans(sJac(:,:,gate2),bin_factor2,3);
m_diff = mB_g-mS_g;
Jac = bJac_g-sJac_g;
m = permute(m_diff, [2 1]);
Jac = permute(Jac, [3 1 2]);
m = m(:);
J = reshape(Jac, [size(Jac,1)*size(Jac,2),size(Jac,3)]);
