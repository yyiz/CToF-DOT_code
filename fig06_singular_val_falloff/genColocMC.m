clear; close all; clc;

% Generate CToF-DOT measurements from single source-detector pair

fdir = "12_12_20_coloc_Jacobian_5x5_singular_val";
fname = "J";

defaultPath = "../dat";

load(sprintf("%s/%s/%s", defaultPath, fdir, fname));

addpath("../lib");

nSD = 25;
shiftPts = linspace(nSD/2, -nSD/2, nSD);
[X, Y] = meshgrid(1:Jheaders.VOX_W, 1:Jheaders.VOX_L);

interpJ = zeros(nSD*nSD*Jheaders.NBINS, nSD*nSD*Jheaders.VOX_H);
for row = 1:nSD
    for col = 1:nSD
        for T = 1:Jheaders.NBINS
            JT = reshape(J(T,:), Jheaders.VOX_L, Jheaders.VOX_W, Jheaders.VOX_H);
            rowShift = shiftPts(row);
            colShift = shiftPts(col);
            
            X_shift = X + colShift;
            Y_shift = Y + rowShift;
            
            JT_interp = zeros(nSD, nSD, Jheaders.VOX_H);
            for z = 1:Jheaders.VOX_H
                indStart = ceil((Jheaders.VOX_L - nSD)/2) + 1;
                indEnd = indStart + nSD - 1;
                JT_interp_z = interp2(X, Y, JT(:,:,z), X_shift, Y_shift);
                JT_interp(:,:,z) = JT_interp_z(indStart:indEnd, indStart:indEnd);
            end
            J_ind = sub2ind([Jheaders.NBINS, nSD, nSD], T, col, row);
            interpJ(J_ind, :) = JT_interp(:)';
        end
    end
    fprintf("Done with row %d/%d\n", row, nSD);
end

figure();
imagesc(interpJ);
set(gca,'ColorScale','log');

S_coloc_25x25 = svd(interpJ);
S_coloc_25x25 = S_coloc_25x25 ./ max(S_coloc_25x25(:));

figure();
semilogy(S_coloc_25x25);

savename = sprintf("%s_coloc", fname);
savepath = sprintf("%s/%s/%s.mat", defaultPath, fdir, savename);
if (exist(savepath, 'file'))
    fprintf("File exists, overwrite?\n");
    keyboard;
end
save(savepath, '-v7.3', 'S_coloc_25x25', 'interpJ', 'J', 'nSD',...
    'shiftPts', 'fdir', 'fname');

