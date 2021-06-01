% Scattering media parameters
U_A = 0.0;
U_S = 0.9;
g = 0;
U_S_red = (1-g)*U_S;
D = 1/(3*(U_A+U_S_red));
c = 0.3; % Speed of light: mm / ps

MFP = 1/(U_A + U_S);
ALB = U_S/(U_A + U_S);

% Timing parameters
compressTime = true;
TMIN = 0;
TMAX = 650;
NBINS = 10;
DTIME = (TMAX - TMIN)/NBINS;
% Perform fluence rate calculation using middle time within each bin
timeAx = transpose((TMIN+(DTIME/2)):DTIME:TMAX);
convL = 2*NBINS-1;

VOX_ORIG = VOX_ORIG - VOXDIM./2; % Set the center of voxel to be VOX_ORIG

paramNames = ["U_A", "U_S", "g", "U_S_red", "D", "coloc",...
    "c", "MFP", "ALB", "SRC_ORIG", "SRC_L", "SRC_W", "SRC_SEP",...
    "DET_ORIG", "DET_L", "DET_W", "DET_SEP", "DET_LEN", "TMIN",...
    "TMAX", "NBINS", "DTIME", "timeAx", "convL",...
    "VOX_ORIG", "VOX_L", "VOX_W", "VOX_H", "VOXDIM"];

for i = 1:length(paramNames)
    params.(paramNames(i)) = eval(paramNames(i));
end

addpath("../lib");

%% Plot scene

if (previewScene)
    plotSetup(params);
    keyboard;
end

%% Generate Jacobian

[J, Jheaders] = genJ(params);

%% Generate background measurements

if (exist('coloc', 'var') == 1 && coloc)
    [sx, sy] = ndgrid(0:SRC_W-1, 0:SRC_L-1);
    [dx, dy] = ndgrid(0:DET_W-1, 0:DET_L-1);
    
    srcPos = zeros(SRC_W, SRC_L, 3);
    srcPos(:,:,1) = SRC_ORIG(1) + sx*SRC_SEP;
    srcPos(:,:,2) = SRC_ORIG(2) + sy*SRC_SEP;
    srcPos(:,:,3) = SRC_ORIG(3);
    srcPos = transpose(reshape(srcPos, prod([SRC_L SRC_W]), []));
    
    r_DET = zeros(SRC_W, SRC_L, 3);
    r_DET(:,:,1) = DET_ORIG(1) + dx*DET_SEP;
    r_DET(:,:,2) = DET_ORIG(2) + dy*DET_SEP;
    r_DET(:,:,3) = DET_ORIG(3);
    r_DET = transpose(reshape(r_DET, prod([DET_L DET_W]), []));
else
    [sx, sy, dx, dy] = ndgrid(0:SRC_W-1, 0:SRC_L-1, 0:DET_W-1, 0:DET_L-1);
    
    srcPos = zeros(SRC_W, SRC_L, DET_W, DET_L, 3);
    srcPos(:,:,:,:,1) = SRC_ORIG(1) + sx*SRC_SEP;
    srcPos(:,:,:,:,2) = SRC_ORIG(2) + sy*SRC_SEP;
    srcPos(:,:,:,:,3) = SRC_ORIG(3);
    srcPos = transpose(reshape(srcPos, prod([SRC_L SRC_W DET_L DET_W]), []));

    r_DET = zeros(SRC_W, SRC_L, DET_W, DET_L, 3);
    r_DET(:,:,:,:,1) = DET_ORIG(1) + dx*DET_SEP;
    r_DET(:,:,:,:,2) = DET_ORIG(2) + dy*DET_SEP;
    r_DET(:,:,:,:,3) = DET_ORIG(3);
    r_DET = transpose(reshape(r_DET, prod([SRC_L SRC_W DET_L DET_W]), []));
end

r_SRC = srcPos; r_SRC(3,:) = r_SRC(3,:) + MFP;
r_IM = srcPos; r_IM(3,:) = r_IM(3,:) - MFP - 4*D;
    
d_sd = vecnorm(r_DET - r_SRC, 2, 1);
d_id = vecnorm(r_DET - r_IM, 2, 1);

Dct = D*c.*timeAx;

refl_sd = (D*c*r_SRC(3, :)./(16*(pi^1.5)*(Dct.^2.5))).*exp(-(d_sd.^2)./(4.*Dct) - U_A*c*timeAx);
refl_id = (D*c*(r_SRC(3, :)+4*D)./(16*(pi^1.5)*(Dct.^2.5))).*exp(-(d_id.^2)./(4.*Dct) - U_A*c*timeAx);

bkgTpsf = (refl_sd + refl_id)*(DET_LEN^2)*DTIME;





