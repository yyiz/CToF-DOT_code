nSrcCol = Jheaders.SRC_W;
nSrcRow = Jheaders.SRC_L;
nDetCol = Jheaders.SENS_W;
nDetRow = Jheaders.SENS_L;
nbins = Jheaders.NBINS;

selInds = zeros(nSrcCol, nSrcRow, nDetCol, nDetRow, nbins);
numMeasPts = prod([nSrcCol, nSrcRow, nDetCol, nDetRow, nbins]);

selInds(srcColInds, srcRowInds, detColInds, detRowInds, :) = 1;
if (selectColoc)
    selInds = zeros(nSrcCol, nSrcRow, nDetCol, nDetRow, nbins);
    for r = srcRowInds
        for c = srcColInds
            selInds(c, r, c, r, :) = 1;
        end
    end
end

indexVec = zeros(numMeasPts, 1);

for sr = 1:nSrcRow
    for sc = 1:nSrcCol
        for dr = 1:nDetRow
            for dc = 1:nDetCol
                
                sc_C = sc-1; sr_C = sr-1; dc_C = dc-1; dr_C = dr-1;
                srcI = sc_C+(sr_C*nSrcCol);
                detI = dc_C+(dr_C*nDetCol);
                datInd = srcI+prod([nSrcCol, nSrcRow])*detI;
                datStart = (datInd * nbins) + 1;
                datEnd = datStart + nbins - 1;
                
                tpsfI = squeeze(selInds(sc,sr,dc,dr,:));
                indexVec(datStart:datEnd) = tpsfI;
            end
        end
    end
end

indexVec = logical(indexVec);

J = J(indexVec, :);

if (~selectColoc)
    Jheaders.SRC_L = length(srcRowInds);
    Jheaders.SRC_W = length(srcColInds);
    Jheaders.SENS_L = length(detRowInds);
    Jheaders.SENS_W = length(detColInds);
else
    Jheaders.SRC_L = length(srcRowInds);
    Jheaders.SRC_W = length(srcColInds);
    Jheaders.SENS_L = 1;
    Jheaders.SENS_W = 1;
end
    
bkgTpsf = bkgTpsf(:);
bkgTpsf = bkgTpsf(indexVec);
bkgTpsf = reshape(bkgTpsf, Jheaders.NBINS, []);


if (selectColoc && ~useMultiplex)
    n = size(bkg_data, 1);
    nbins = size(bkg_data, 5);
    nmulti = size(bkg_data, 6);
    ntacq = size(bkg_data, 7);
    temp_meas = zeros(n, n, 1, 1, nbins, nmulti, ntacq);
    temp_bkg = zeros(n, n, 1, 1, nbins, nmulti, ntacq);
    for r = 1:4
        for c = 1:4
            temp_meas(c, r, 1, 1, :) = measure_data(c, r, c, r, :);
            temp_bkg(c, r, 1, 1, :) = bkg_data(c, r, c, r, :);
        end
    end
    measure_data = temp_meas;
    bkg_data = temp_bkg;
end

if size(bkg_data,1) ~= length(srcColInds)
    measure_data = measure_data(srcColInds, srcRowInds, detColInds, detRowInds, :);
    bkg_data = bkg_data(srcColInds, srcRowInds, detColInds, detRowInds, :);
end
