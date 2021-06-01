function [calibM_mc] = calibExpData(expDatFname, bkgM_mc, homo, target,...
                                    s, filts)

load(expDatFname);

nSrcCol = size(homo, 1);
nSrcRow = size(homo, 2);
nDetCol = size(homo, 3);
nDetRow = size(homo, 4);
nbins = size(homo, 5);

expTimeAx=1:nbins;

n_filts = length(filts);
n_sd = nSrcCol*nSrcRow*nDetCol*nDetRow;
sizeDat = n_sd * n_filts;
calibM_mc = zeros(sizeDat, 1);

useEfilt = contains(filts, 'E');
useMfilt = contains(filts, 'M');
useLfilt = contains(filts, 'L');

bkgE_pos = strfind(filts, 'E') - 1;
bkgM_pos = strfind(filts, 'M') - 1;
bkgL_pos = strfind(filts, 'L') - 1;

for sr = 1:nSrcRow
    for sc = 1:nSrcCol
        for dr = 1:nDetRow
            for dc = 1:nDetCol

                homo_raw = squeeze(homo(sc,sr,dc,dr,:));
                tar_raw = squeeze(target(sc,sr,dc,dr,:));
    
                smoothed_homo=smoothdata(homo_raw,'gaussian',100);
                smoothed_tar=smoothdata(tar_raw,'gaussian',100);
                
                sc_C = sc-1; sr_C = sr-1; dc_C = dc-1; dr_C = dr-1;
                srcI = sc_C+(sr_C*nSrcCol);
                detI = dc_C+(dr_C*nDetCol);
                bkgInd = srcI+prod([nSrcCol, nSrcRow])*detI + 1;
                
                if (useEfilt)
                    homo_E = trapz(expTimeAx,smoothed_homo');        
                    tar_E = trapz(expTimeAx,smoothed_tar');
                    bkg_E_ind = (bkgE_pos * n_sd) + bkgInd;
                    calibM_mc(bkg_E_ind) = tar_E*bkgM_mc(bkg_E_ind)/homo_E;
                end
                if (useMfilt)
                    homo_M = trapz(expTimeAx,(expTimeAx.*smoothed_homo'));
                    tar_M = trapz(expTimeAx,(expTimeAx.*smoothed_tar'));
                    bkg_M_ind = (bkgM_pos * n_sd) + bkgInd;
                    calibM_mc(bkg_M_ind) = tar_M*bkgM_mc(bkg_M_ind)/homo_M;
                end
                if (useLfilt)
                    homo_L = trapz(expTimeAx,(exp(-s*expTimeAx).*smoothed_homo'));
                    tar_L = trapz(expTimeAx,(exp(-s*expTimeAx).*smoothed_tar'));
                    bkg_L_ind = (bkgL_pos * n_sd) + bkgInd;
                    calibM_mc(bkg_L_ind) = tar_L*bkgM_mc(bkg_L_ind)/homo_L;
                end

            end
        end
    end
end

end