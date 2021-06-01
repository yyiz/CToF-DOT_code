function [J, Jheaders] = genJ(params)

loadParams(params);

%% Determine Jacobian
if (coloc)
    J = zeros(SRC_L*SRC_W*NBINS, VOX_L*VOX_W*VOX_H);
else
    J = zeros(SRC_L*SRC_W*DET_L*DET_W*NBINS, VOX_W*VOX_L*VOX_H);
end

[voxXi, voxYi, voxZi] = ndgrid(0:VOX_W-1, 0:VOX_L-1, 0:VOX_H-1);

voxCoords = zeros(VOX_W, VOX_L, VOX_H, 3);
voxCoords(:,:,:,1) = VOX_ORIG(1) + VOXDIM(1)*(0.5 + voxXi); % x-coordinates
voxCoords(:,:,:,2) = VOX_ORIG(2) + VOXDIM(2)*(0.5 + voxYi); % y-coordinates
voxCoords(:,:,:,3) = VOX_ORIG(3) + VOXDIM(3)*(0.5 + voxZi); % z-coordinates

VOX = transpose(reshape(voxCoords, prod([VOX_L, VOX_W, VOX_H]), []));

if (coloc)
    drRange = 1;
    dcRange = 1;
else
    drRange = 1:DET_L;
    dcRange = 1:DET_W;
end

tic;
for sr = 1:SRC_L
for sc = 1:SRC_W
for dr = drRange
for dc = dcRange

r_SRC = SRC_ORIG + [(sc-1)*SRC_SEP; (sr-1)*SRC_SEP; 0];
if (coloc)
    r_DET = DET_ORIG + [(sc-1)*SRC_SEP; (sr-1)*SRC_SEP; 0];
else
    r_DET = DET_ORIG + [(dc-1)*SRC_SEP; (dr-1)*SRC_SEP; 0];
end
    
r_MFP = [r_SRC(1); r_SRC(2); MFP];
r_IMS = [r_SRC(1); r_SRC(2); -MFP - 4*D];
r_IMD = [VOX(1,:); VOX(2,:); -VOX(3,:)-4*D];

%% Source path estimate

d_SRC_VOX = vecnorm(r_MFP - VOX, 2, 1);
d_IMS_VOX = vecnorm(r_IMS - VOX, 2, 1);

green_SRC_VOX = (ALB*c./((4*pi*D*c*timeAx).^1.5)).*exp(-(d_SRC_VOX.^2)./(4*D*c.*timeAx) - U_A*c.*timeAx);
green_IMS_VOX = (ALB*c./((4*pi*D*c*timeAx).^1.5)).*exp(-(d_IMS_VOX.^2)./(4*D*c.*timeAx) - U_A*c.*timeAx);

green_SRC_VOX = (green_SRC_VOX - green_IMS_VOX) * DTIME * prod(VOXDIM);

%% Detector path estimate

d_VOX_DET = vecnorm(r_DET - VOX, 2, 1);
d_IMD_DET = vecnorm(r_DET - r_IMD, 2, 1);

green_VOX_DET = (D*c*VOX(3,:)./(16*(pi^1.5)*((D*c.*timeAx).^2.5))).*exp(-(d_VOX_DET.^2)./(4.*(D*c.*timeAx)) - U_A*c*timeAx);
green_IMD_DET = (D*c*(VOX(3,:)+4*D)./(16*(pi^1.5)*((D*c.*timeAx).^2.5))).*exp(-(d_IMD_DET.^2)./(4.*(D*c.*timeAx)) - U_A*c*timeAx);

green_VOX_DET = (green_VOX_DET + green_IMD_DET)*(DET_LEN^2)*DTIME;
%% Generate full path: closed form

% implement convolution using fft
tpsf_full = ifft(fft(green_SRC_VOX, convL, 1) .* fft(green_VOX_DET, convL, 1), convL, 1);

if (coloc)
    detInd = 0;
else
    detInd = (dr-1)*DET_W + (dc-1);
end
srcInd = (sr-1)*SRC_W + (sc-1);
srcDetInd = (detInd * (SRC_L*SRC_W)) + srcInd;
Jind1 = srcDetInd*NBINS + 1;
Jind2 = Jind1 + NBINS - 1;

J(Jind1:Jind2,:) = tpsf_full(1:NBINS,:);

end
end
end
end
genJTime = toc;
fprintf("Generate Jacobian Runtime: %.3f seconds\n", genJTime);
%% Save variables

Jheaders.SRC_L = SRC_L;
Jheaders.SRC_W = SRC_W;
Jheaders.SENS_L = DET_L;
Jheaders.SENS_W = DET_W;
Jheaders.VOX_L = VOX_L;
Jheaders.VOX_W = VOX_W;
Jheaders.VOX_H = VOX_H;
Jheaders.VOX_ORIGX = VOX_ORIG(1);
Jheaders.VOX_ORIGY = VOX_ORIG(2);
Jheaders.VOX_ORIGZ = VOX_ORIG(3);
Jheaders.VOX_SIDELEN = VOXDIM(1);
Jheaders.VOX_ZLEN = VOXDIM(3);
Jheaders.NBINS = NBINS;
Jheaders.TIME_MIN = TMIN;
Jheaders.TIME_MAX = TMAX;

end

function loadParams(params)
    vars = fieldnames(params);
    for i = 1:length(vars)
        assignin('caller', vars{i}, params.(vars{i}));
    end
end
