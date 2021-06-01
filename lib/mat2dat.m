function matFull = mat2dat(matIn, Jheaders, isJ)

nSrcCol = Jheaders.SRC_W;
nSrcRow = Jheaders.SRC_L;
nDetCol = Jheaders.SENS_W;
nDetRow = Jheaders.SENS_L;
nbins = Jheaders.NBINS;
nvox = prod([Jheaders.VOX_L, Jheaders.VOX_W, Jheaders.VOX_H]);

if isJ
    
    matFull = zeros(nSrcCol, nSrcRow, nDetCol, nDetRow, nbins, nvox);
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

                    matFull(sc, sr, dc, dr,:,:) = matIn(datStart:datEnd,:);
                end
            end
        end
    end
    
else

    matFull = zeros(nSrcCol, nSrcRow, nDetCol, nDetRow, nbins);
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

                    matFull(sc, sr, dc, dr,:) = matIn(datStart:datEnd);
                end
            end
        end
    end

end

end