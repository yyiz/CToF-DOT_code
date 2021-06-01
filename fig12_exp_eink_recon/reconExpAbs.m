%% Set up reconstruction parameters

filts = "M";
nFilts = strlength(filts);
s = 0.0001;

useMultiplex = false;
cropJacobian = true;
selectColoc =  false;
srcRowInds = 1:10;
srcColInds = 1:10;
detRowInds = 1:10;
detColInds = 1:10;

gateExpDat = false;
expTimeGate = 5700:9000;

repInd = 1;
waveletLvl = 3;
domain = "mvp";

defaultPath = "../dat";
Jfname = sprintf("%s/9_21_20_jacobian_10x10_skull_6p5/J.mat", defaultPath);
bkgFname = sprintf("%s/9_21_20_jacobian_10x10_skull_6p5/tpsf0.mat", defaultPath);
expDatFname = sprintf("%s/02-Dec-2020_eink_meas/%s.mat", defaultPath, expFname);
truthFile = sprintf("%s/truth_ims/%s", defaultPath, truthImFile);

showTruth = true;

shouldSave = true;
savepath = "../../figs";
saveDir = "fig12_exp_eink_recon";
savePath = sprintf("%s/%s", savepath, saveDir);
saveName = expFname;

%% Pre-process data

addpath("../lib/");

% Load (simulated) Jacobian and background TPSF
load(Jfname);
load(bkgFname);
load(expDatFname);

if (size(bkg_data, 6) > 1)
    bkg_data = bkg_data(:,:,:,:,:,repInd);
    measure_data = measure_data(:,:,:,:,:,repInd);
end

if (cropJacobian)
    selectSrcDet;
end
[J_mc, bkgM_mc] = tempFilts(J, bkgTpsf, Jheaders, filts, 's', s);

% Calibrate data to bkg 
if (gateExpDat)
    bkg_data = bkg_data(:,:,:,:,expTimeGate);
    measure_data = measure_data(:,:,:,:,expTimeGate);
end
[calibM_mc] = calibExpData(expDatFname, bkgM_mc, bkg_data, measure_data,...
                           s, filts);
m = bkgM_mc - calibM_mc;


if (strcmp(domain, "mvp"))
    sizeX = size(J_mc,2);
elseif (endsWith(domain, "dct"))
    sizeX = [Jheaders.VOX_L, Jheaders.VOX_W];
    if (fistaOpts.nonneg)
        warning("Set fistaOpts.nonneg=false when using DCT transform");
        keyboard;
    end
elseif (endsWith(domain, "dwt"))
    sizeX = [Jheaders.VOX_L, Jheaders.VOX_W];
    if (fistaOpts.nonneg)
        warning("Set fistaOpts.nonneg=false when using Wavelet transform");
        keyboard;
    end
else
    error("invalid transform domain type");
end
[f_fista, fT_fista, fistaStep] = mat2Handle(J_mc, domain,...
    'VOX_L', Jheaders.VOX_L, 'VOX_W', Jheaders.VOX_W, 'lvl', waveletLvl);


%% Perform image reconstruction

fprintf("Running FISTA...");
tic; [fistaRecon, fistaErr] = fista(m,f_fista,fT_fista,fistaStep,sizeX,fistaOpts); fistaRuntime = toc;
fprintf("done! Finished in: %.3f seconds\n", fistaRuntime);


%% Post-processing and plotting

% Plot FISTA image reconstruction
fistaReconIm = reshape(fistaRecon, [Jheaders.VOX_L, Jheaders.VOX_W, Jheaders.VOX_H]);
fistaReconIm = rot90(fistaReconIm, 1);

f1 = figure();
truthIm = imread(truthFile);
imagesc(truthIm);
axis image;
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

f2 = figure();
imagesc(fistaReconIm);
axis image;
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

%% Save results

if ~(exist(savePath, 'dir'))
    mkdir(savePath);
end

export_fig(f1, sprintf("%s/reconExpAbs_%s_truth.png", savePath, saveName), '-m3', '-transparent', '-png');
export_fig(f2, sprintf("%s/reconExpAbs_%s_recon.png", savePath, saveName), '-m3', '-transparent', '-png');




