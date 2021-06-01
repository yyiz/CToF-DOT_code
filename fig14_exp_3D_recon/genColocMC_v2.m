%% Calculate full collocated Jacobian matrix

shiftPts_L = linspace(sdArr_Len/2, -sdArr_Len/2, nSD_L);
shiftPts_W = linspace(sdArr_Width/2, -sdArr_Width/2, nSD_W);

[X, Y] = meshgrid(1:VOX_W, 1:VOX_L);

J_final = zeros(nSD_L*nSD_W*NBINS, voxSize_L*voxSize_W*VOX_H);
for row = 1:nSD_L
    for col = 1:nSD_W
        for T = 1:NBINS
            JT = reshape(J_postproc(T,:), VOX_L, VOX_W, VOX_H);
            rowShift = shiftPts_L(row);
            colShift = shiftPts_W(col);
            
            X_shift = X + colShift;
            Y_shift = Y + rowShift;
            
            JT_interp = zeros(voxSize_L, voxSize_W, VOX_H);
            for z = 1:VOX_H
                indStart_1 = ceil((VOX_L - voxSize_L)/2) + 1;
                indEnd_1 = indStart_1 + voxSize_L - 1;
                
                indStart_2 = ceil((VOX_W - voxSize_W)/2) + 1;
                indEnd_2 = indStart_2 + voxSize_W - 1;
                
                JT_interp_z = interp2(X, Y, JT(:,:,z), X_shift, Y_shift, 'spline');
                JT_interp(:,:,z) = JT_interp_z(indStart_1:indEnd_1, indStart_2:indEnd_2);
            end
            J_ind = sub2ind([NBINS, nSD_L, nSD_W], T, row, col);      
            
            JT_interp(JT_interp < 0) = 0;
            JT_interp = JT_interp .^ gamma;
            
            J_final(J_ind, :) = JT_interp(:);
            
        end
    end
    fprintf("Done with row %d/%d\n", row, nSD_L);
end

SRC_L = nSD_L;
SRC_W = nSD_W;
DET_L = 1; % Because colocated
DET_W = 1;
VOX_L = voxSize_L;
VOX_W = voxSize_W;

Jheaders.SRC_L = nSD_L;
Jheaders.SRC_W = nSD_W;
Jheaders.DET_L = 1; % Because colocated
Jheaders.DET_W = 1;
Jheaders.VOX_L = voxSize_L;
Jheaders.VOX_W = voxSize_W;


